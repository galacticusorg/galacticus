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

module fodeiv2
  use FGSL
  use, intrinsic :: iso_c_binding
  implicit none
  private
! ordinary differential equations
  public :: fodeiv2_system_init, fodeiv2_system_free, &
       fodeiv2_step_alloc, fodeiv2_step_status, fodeiv2_system_status, &
       fodeiv2_step_reset, fodeiv2_step_free, fodeiv2_step_name, &
       fodeiv2_step_order, fodeiv2_step_apply, fodeiv2_control_standard_new, &
       fodeiv2_control_y_new, fodeiv2_control_yp_new, &
       fodeiv2_control_scaled_new, fodeiv2_control_init, &
       fodeiv2_control_free, fodeiv2_control_hadjust, &
       fodeiv2_control_name, fodeiv2_evolve_alloc, fodeiv2_evolve_apply, &
       fodeiv2_evolve_reset, fodeiv2_evolve_free, &
       fodeiv2_driver_apply, fodeiv2_driver_reset, &
       fodeiv2_driver_free, fodeiv2_driver_alloc_standard_new, &
       fodeiv2_driver_alloc_y_new, fodeiv2_driver_alloc_yp_new, &
       fodeiv2_driver_alloc_scaled_new, fodeiv2_driver_set_hmin, &
       fodeiv2_driver_set_hmax, fodeiv2_driver_set_nmax, fodeiv2_driver_status

  ! Specify an explicit dependence on the bivar.o object file.
  !: ./work/build/numerical.ODE_solver.ODEIV2.utils.o

!
! Types: Ordinary Differential Equations
!
  type, public :: fodeiv2_system
     private
     type(c_ptr) :: odeiv2_system=c_null_ptr 
  end type fodeiv2_system
  type, public :: fodeiv2_step_type
     private
     integer(kind=c_int) :: which=0 
  end type fodeiv2_step_type
  type  (fodeiv2_step_type), parameter, public :: fodeiv2_step_bsimp=fodeiv2_step_type(8) , fodeiv2_step_msadams=fodeiv2_step_type(10), & 
       &                                          fodeiv2_step_msbdf=fodeiv2_step_type(11), fodeiv2_step_rk1imp =fodeiv2_step_type(9) , & 
       &                                          fodeiv2_step_rk2  =fodeiv2_step_type(1) , fodeiv2_step_rk2imp =fodeiv2_step_type(6) , & 
       &                                          fodeiv2_step_rk4  =fodeiv2_step_type(2) , fodeiv2_step_rk4imp =fodeiv2_step_type(7) , & 
       &                                          fodeiv2_step_rk8pd=fodeiv2_step_type(5) , fodeiv2_step_rkck   =fodeiv2_step_type(4) , & 
       &                                          fodeiv2_step_rkf45=fodeiv2_step_type(3)                                                 
  
  type, public :: fodeiv2_step
     type(c_ptr) :: odeiv2_step=c_null_ptr 
  end type fodeiv2_step
  type, public :: fodeiv2_control
     type(c_ptr) :: odeiv2_control=c_null_ptr 
  end type fodeiv2_control
  integer(kind=fgsl_int), parameter, public :: fodeiv2_hadj_inc=1  
  integer(kind=fgsl_int), parameter, public :: fodeiv2_hadj_nil=0  
  integer(kind=fgsl_int), parameter, public :: fodeiv2_hadj_dec=-1 
  type, public :: fodeiv2_evolve
     type(c_ptr) :: odeiv2_evolve 
  end type fodeiv2_evolve
  type, public :: fodeiv2_driver
     type(c_ptr) :: odeiv2_driver 
  end type fodeiv2_driver

! required C interfaces
! FGSL names occurring here are auxiliary routines
! needed to transfer static C information to the Fortran subsystem
  interface
!-*-f90-*-
!
!  Interfaces: Ordinary differential equations
!
     function gsl_odeiv2_step_alloc(t, dim) bind(c)
       import
       type   (c_ptr        ), value :: t                     
       integer(kind=c_size_t), value :: dim                   
       type   (c_ptr        )        :: gsl_odeiv2_step_alloc 
     end function gsl_odeiv2_step_alloc
     function gsl_odeiv2_aux_odeiv_step_alloc(step_type) bind(c)
       import
       integer(kind=c_int), value :: step_type                       
       type   (c_ptr     )        :: gsl_odeiv2_aux_odeiv_step_alloc 
     end function gsl_odeiv2_aux_odeiv_step_alloc
     function fodeiv2_system_cinit(func, dimension, params, jacobian) bind(c)
       import
       type   (c_funptr     ), value :: func                 
       integer(kind=c_size_t), value :: dimension            
       type   (c_ptr        ), value :: params               
       type   (c_funptr     ), value :: jacobian             
       type   (c_ptr        )        :: fodeiv2_system_cinit 
     end function fodeiv2_system_cinit
     subroutine fodeiv2_system_cfree(system) bind(c)
       import
       type(c_ptr), value :: system 
     end subroutine fodeiv2_system_cfree
     function gsl_odeiv2_step_reset(s) bind(c)
       import 
       type   (c_ptr     ), value :: s                     
       integer(kind=c_int)        :: gsl_odeiv2_step_reset 
     end function gsl_odeiv2_step_reset
     subroutine gsl_odeiv2_step_free(s) bind(c)
       import 
       type(c_ptr), value :: s 
     end subroutine gsl_odeiv2_step_free
     function gsl_odeiv2_step_name (s) bind(c)
       import
       type(c_ptr), value :: s                    
       type(c_ptr)        :: gsl_odeiv2_step_name 
     end function gsl_odeiv2_step_name
     function gsl_odeiv2_step_order(s) bind(c)
       import
       type   (c_ptr     ), value :: s                     
       integer(kind=c_int)        :: gsl_odeiv2_step_order 
     end function gsl_odeiv2_step_order
     function gsl_odeiv2_step_apply(s, t, h, y, yerr, dydt_in, dydt_out, dydt) bind(c)
       import
       type   (c_ptr        ), value                       :: s                                     
       real   (kind=c_double), value                       :: h                    , t              
       real   (kind=c_double), dimension(*), intent(inout) :: dydt_in              , dydt_out, y, & 
            &                                                 yerr                                  
       type   (c_ptr        ), value                       :: dydt                                  
       integer(kind=c_int   )                              :: gsl_odeiv2_step_apply                 
     end function gsl_odeiv2_step_apply
     function gsl_odeiv2_control_standard_new(eps_abs, eps_rel, a_y, a_dydt) bind(c)
       import
       real  (kind=c_double), value :: a_dydt                         , a_y    , & 
            &                          eps_abs                        , eps_rel    
       type  (c_ptr        )        :: gsl_odeiv2_control_standard_new             
     end function gsl_odeiv2_control_standard_new
     function gsl_odeiv2_control_y_new(eps_abs, eps_rel) bind(c)
       import
       real(kind=c_double), value :: eps_abs                 , eps_rel 
       type(c_ptr        )        :: gsl_odeiv2_control_y_new          
     end function gsl_odeiv2_control_y_new
     function gsl_odeiv2_control_yp_new(eps_abs, eps_rel) bind(c)
       import
       real(kind=c_double), value :: eps_abs                  , eps_rel 
       type(c_ptr        )        :: gsl_odeiv2_control_yp_new          
     end function gsl_odeiv2_control_yp_new
     function gsl_odeiv2_control_scaled_new(eps_abs, eps_rel, a_y, a_dydt, &
          scale_abs, dim) bind(c)
       import
       real   (kind=c_double), value                       :: a_dydt                       , a_y    , & 
            &                                                 eps_abs                      , eps_rel    
       real   (kind=c_double), dimension(*), intent(in   ) :: scale_abs                                 
       integer(kind=c_size_t), value                       :: dim                                       
       type   (c_ptr        )                              :: gsl_odeiv2_control_scaled_new             
     end function gsl_odeiv2_control_scaled_new
! odeiv2_control_alloc presently not attached
     function gsl_odeiv2_control_init(c, eps_abs, eps_rel, a_y, a_dydt) bind(c)
       import
       type   (c_ptr        ), value :: c                                   
       real   (kind=c_double), value :: a_dydt                 , a_y    , & 
            &                           eps_abs                , eps_rel    
       integer(kind=c_int   )        :: gsl_odeiv2_control_init             
     end function gsl_odeiv2_control_init
     subroutine gsl_odeiv2_control_free(c) bind(c)
       import
       type(c_ptr), value :: c 
     end subroutine gsl_odeiv2_control_free
     function gsl_odeiv2_control_hadjust(c, s, y0, yerr, dydt, h) bind(c)
       import
       type   (c_ptr        ), value                       :: c                         , s     
       real   (kind=c_double), dimension(*), intent(in   ) :: dydt                      , y0, & 
            &                                                 yerr                              
       real   (kind=c_double), dimension(*), intent(inout) :: h                                 
       integer(kind=c_int   )                              :: gsl_odeiv2_control_hadjust        
     end function gsl_odeiv2_control_hadjust
     function gsl_odeiv2_control_name (s) bind(c)
       import
       type(c_ptr), value :: s                       
       type(c_ptr)        :: gsl_odeiv2_control_name 
     end function gsl_odeiv2_control_name
     function gsl_odeiv2_evolve_alloc(dim) bind(c)
       import
       integer(kind=c_size_t), value :: dim                     
       type   (c_ptr        )        :: gsl_odeiv2_evolve_alloc 
     end function gsl_odeiv2_evolve_alloc
     function gsl_odeiv2_evolve_apply(e, con, step, dydt, t, t1, h, y) bind(c)
       import
       type   (c_ptr        ), value                       :: con                    , dydt, & 
            &                                                 e                      , step    
       real   (kind=c_double), dimension(*), intent(inout) :: y                                
       real   (kind=c_double)              , intent(inout) :: h                      , t       
       real   (kind=c_double), value                       :: t1                               
       integer(kind=c_int   )                              :: gsl_odeiv2_evolve_apply          
     end function gsl_odeiv2_evolve_apply
     function gsl_odeiv2_evolve_reset(s) bind(c)
       import 
       type   (c_ptr     ), value :: s                       
       integer(kind=c_int)        :: gsl_odeiv2_evolve_reset 
     end function gsl_odeiv2_evolve_reset
     subroutine gsl_odeiv2_evolve_free(s) bind(c)
       import 
       type(c_ptr), value :: s 
     end subroutine gsl_odeiv2_evolve_free
     function gsl_odeiv2_driver_apply(d, t, t1, y &
#ifdef PROFILE
          &                           , sa        &
#endif
          &                          ) bind(c)
       import
       type(c_ptr        ), value                       :: d  
       real(kind=c_double)                              :: t  
       real(kind=c_double), value                       :: t1 
       real(kind=c_double), dimension(*), intent(inout) :: y  
#ifdef PROFILE
       type(c_funptr), value :: sa 
#endif
       integer(kind=c_int) :: gsl_odeiv2_driver_apply 
     end function gsl_odeiv2_driver_apply
     function gsl_odeiv2_driver_reset(d) bind(c)
       import 
       type   (c_ptr     ), value :: d                       
       integer(kind=c_int)        :: gsl_odeiv2_driver_reset 
     end function gsl_odeiv2_driver_reset
     subroutine gsl_odeiv2_driver_free(d) bind(c)
       import 
       type(c_ptr), value :: d 
     end subroutine gsl_odeiv2_driver_free
     function gsl_odeiv2_driver_alloc_standard_new(sys, t, hstart, eps_abs, eps_rel, a_y, a_dydt) bind(c)
       import
       type  (c_ptr        ), value :: sys                                 , t          
       real  (kind=c_double), value :: a_dydt                              , a_y    , & 
            &                          eps_abs                             , eps_rel, & 
            &                          hstart                                           
       type  (c_ptr        )        :: gsl_odeiv2_driver_alloc_standard_new             
     end function gsl_odeiv2_driver_alloc_standard_new
     function gsl_odeiv2_driver_alloc_y_new(sys, t, hstart, eps_abs, eps_rel) bind(c)
       import
       type  (c_ptr        ), value :: sys                          , t          
       real  (kind=c_double), value :: eps_abs                      , eps_rel, & 
            &                          hstart                                    
       type  (c_ptr        )        :: gsl_odeiv2_driver_alloc_y_new             
     end function gsl_odeiv2_driver_alloc_y_new
     function gsl_odeiv2_driver_alloc_yp_new(sys, t, hstart, eps_abs, eps_rel) bind(c)
       import
       type  (c_ptr        ), value :: sys                           , t          
       real  (kind=c_double), value :: eps_abs                       , eps_rel, & 
            &                          hstart                                     
       type  (c_ptr        )        :: gsl_odeiv2_driver_alloc_yp_new             
     end function gsl_odeiv2_driver_alloc_yp_new
     function gsl_odeiv2_driver_alloc_scaled_new(sys, t, hstart, eps_abs, eps_rel, a_y, a_dydt, &
          scale_abs) bind(c)
       import
       type  (c_ptr        ), value                       :: sys                               , t          
       real  (kind=c_double), value                       :: a_dydt                            , a_y    , & 
            &                                                eps_abs                           , eps_rel, & 
            &                                                hstart                                         
       real  (kind=c_double), dimension(*), intent(in   ) :: scale_abs                                      
       type  (c_ptr        )                              :: gsl_odeiv2_driver_alloc_scaled_new             
     end function gsl_odeiv2_driver_alloc_scaled_new
     function gsl_odeiv2_driver_set_hmin(d, hmin) bind(c)
       import
       type   (c_ptr        ), value :: d                          
       real   (kind=c_double), value :: hmin                       
       integer(kind=c_int   )        :: gsl_odeiv2_driver_set_hmin 
     end function gsl_odeiv2_driver_set_hmin
     function gsl_odeiv2_driver_set_hmax(d, hmax) bind(c)
       import
       type   (c_ptr        ), value :: d                          
       real   (kind=c_double), value :: hmax                       
       integer(kind=c_int   )        :: gsl_odeiv2_driver_set_hmax 
     end function gsl_odeiv2_driver_set_hmax
     function gsl_odeiv2_driver_set_nmax(d, nmax) bind(c)
       import
       type   (c_ptr      ), value :: d                          
       integer(kind=c_long), value :: nmax                       
       integer(kind=c_int )        :: gsl_odeiv2_driver_set_nmax 
     end function gsl_odeiv2_driver_set_nmax
  end interface
contains
!-*-f90-*-
!
! API: Ordinary differential equations
!
  function fodeiv2_system_init(func, dimension, params, jacobian)
    optional :: jacobian
    interface
       function func(t, y, dydt, params) bind(c)
         use, intrinsic :: iso_c_binding
         real   (kind=c_double)              , value         :: t      
         real   (kind=c_double), dimension(*), intent(in   ) :: y      
         real   (kind=c_double), dimension(*)                :: dydt   
         type   (c_ptr        )              , value         :: params 
         integer(kind=c_int   )                              :: func   
       end function func
       function jacobian(t, y, dfdy, dfdt, params) bind(c)
         use, intrinsic :: iso_c_binding
         real   (kind=c_double)              , value         :: t        
         real   (kind=c_double), dimension(*), intent(in   ) :: y        
         real   (kind=c_double), dimension(*)                :: dfdy     
         real   (kind=c_double), dimension(*)                :: dfdt     
         type   (c_ptr        )              , value         :: params   
         integer(kind=c_int   )                              :: jacobian 
       end function jacobian
    end interface
    integer(kind=fgsl_size_t)                          :: dimension           
    type   (c_ptr           ), intent(in   ), optional :: params              
    type   (fodeiv2_system  )                          :: fodeiv2_system_init 
    type   (c_funptr        )                          :: func_loc            
    type   (c_funptr        )                          :: jacobian_loc        
    type   (c_ptr           )                          :: params_loc          
    func_loc = c_funloc(func)
    params_loc = c_null_ptr
    jacobian_loc = c_null_funptr
    if (present(jacobian)) jacobian_loc = c_funloc(jacobian)
    if (present(params)) params_loc = params
    fodeiv2_system_init%odeiv2_system = &
         fodeiv2_system_cinit(func_loc, dimension, &
         params_loc, jacobian_loc)
  end function fodeiv2_system_init
  subroutine fodeiv2_system_free(system)
    type(fodeiv2_system), intent(inout) :: system 
    call fodeiv2_system_cfree(system%odeiv2_system)
  end subroutine fodeiv2_system_free
  function fodeiv2_step_alloc(t, dim) 
    type   (fodeiv2_step_type), intent(in   ) :: t                  
    integer(kind=fgsl_size_t ), intent(in   ) :: dim                
    type   (fodeiv2_step     )                :: fodeiv2_step_alloc 
    !
    type   (c_ptr            )                :: step_type          
    step_type = gsl_odeiv2_aux_odeiv_step_alloc(t%which)
    if (c_associated(step_type)) then
       fodeiv2_step_alloc%odeiv2_step = gsl_odeiv2_step_alloc(step_type, dim)
    else
       fodeiv2_step_alloc%odeiv2_step = c_null_ptr
    end if
  end function fodeiv2_step_alloc
  function fodeiv2_step_reset(s)
    type   (fodeiv2_step ), intent(inout) :: s                  
    integer(kind=fgsl_int)                :: fodeiv2_step_reset 
    fodeiv2_step_reset = gsl_odeiv2_step_reset(s%odeiv2_step)
  end function fodeiv2_step_reset
  subroutine fodeiv2_step_free(s)
    type(fodeiv2_step), intent(inout) :: s 
    call gsl_odeiv2_step_free(s%odeiv2_step)
  end subroutine fodeiv2_step_free
  function fodeiv2_step_name (s)
    type     (fodeiv2_step                   ), intent(in   ) :: s                 
    character(kind=fgsl_char, len=fgsl_strmax)                :: fodeiv2_step_name 
    !
    type     (c_ptr                          )                :: name              
    name = gsl_odeiv2_step_name(s%odeiv2_step)
    fodeiv2_step_name = fgsl_name(name)
  end function fodeiv2_step_name
  function fodeiv2_step_order(s)
    type   (fodeiv2_step ), intent(in   ) :: s                  
    integer(kind=fgsl_int)                :: fodeiv2_step_order 
    fodeiv2_step_order = gsl_odeiv2_step_order(s%odeiv2_step)
  end function fodeiv2_step_order
  function fodeiv2_step_apply(s, t, h, y, yerr, dydt_in, dydt_out, dydt)
    type   (fodeiv2_step    ), intent(in   ) :: s                                           
    real   (kind=fgsl_double), intent(in   ) :: h                    , t                    
    real   (kind=fgsl_double), intent(inout) :: dydt_in           (:), dydt_out(:), y(:), & 
         &                                      yerr              (:)                       
    type   (fodeiv2_system  ), intent(in   ) :: dydt                                        
    integer(kind=fgsl_int   )                :: fodeiv2_step_apply                          
    fodeiv2_step_apply = gsl_odeiv2_step_apply(s%odeiv2_step, t, h, y, yerr, &
         dydt_in, dydt_out, dydt%odeiv2_system)
  end function fodeiv2_step_apply
  function fodeiv2_control_standard_new(eps_abs, eps_rel, a_y, a_dydt)
    real  (kind=fgsl_double), intent(in   ) :: a_dydt                      , a_y    , & 
         &                                     eps_abs                     , eps_rel    
    type  (fodeiv2_control )                :: fodeiv2_control_standard_new             
    fodeiv2_control_standard_new%odeiv2_control = &
         gsl_odeiv2_control_standard_new(eps_abs, eps_rel, a_y, a_dydt)
  end function fodeiv2_control_standard_new
  function fodeiv2_control_y_new(eps_abs, eps_rel)
    real(kind=fgsl_double), intent(in   ) :: eps_abs              , eps_rel 
    type(fodeiv2_control )                :: fodeiv2_control_y_new          
    fodeiv2_control_y_new%odeiv2_control = &
         gsl_odeiv2_control_y_new(eps_abs, eps_rel)
  end function fodeiv2_control_y_new
  function fodeiv2_control_yp_new(eps_abs, eps_rel)
    real(kind=fgsl_double), intent(in   ) :: eps_abs               , eps_rel 
    type(fodeiv2_control )                :: fodeiv2_control_yp_new          
    fodeiv2_control_yp_new%odeiv2_control = &
         gsl_odeiv2_control_yp_new(eps_abs, eps_rel)
  end function fodeiv2_control_yp_new
  function fodeiv2_control_scaled_new(eps_abs, eps_rel, a_y, a_dydt, scale_abs, dim)
    real   (kind=fgsl_double), intent(in   ) :: a_dydt                       , a_y    , & 
         &                                      eps_abs                      , eps_rel    
    real   (kind=fgsl_double), intent(in   ) :: scale_abs                 (:)             
    integer(kind=fgsl_size_t), intent(in   ) :: dim                                       
    type   (fodeiv2_control )                :: fodeiv2_control_scaled_new                
    fodeiv2_control_scaled_new%odeiv2_control = &
         gsl_odeiv2_control_scaled_new(eps_abs, eps_rel, a_y, a_dydt, scale_abs, dim)
  end function fodeiv2_control_scaled_new
! FIXME (?) fodeiv2_control_alloc is presently not implemented
  function fodeiv2_control_init(c, eps_abs, eps_rel, a_y, a_dydt)
    type   (fodeiv2_control ), intent(in   ) :: c                                     
    real   (kind=fgsl_double), intent(in   ) :: a_dydt              , a_y, eps_abs, & 
         &                                      eps_rel                               
    integer(kind=fgsl_int   )                :: fodeiv2_control_init                  
    fodeiv2_control_init = &
         gsl_odeiv2_control_init(c%odeiv2_control, eps_abs, eps_rel, a_y, a_dydt)
  end function fodeiv2_control_init
  subroutine fodeiv2_control_free(c)
    type(fodeiv2_control), intent(inout) :: c 
    call gsl_odeiv2_control_free(c%odeiv2_control)
  end subroutine fodeiv2_control_free
  function fodeiv2_control_hadjust(c, s, y0, yerr, dydt, h) 
    type   (fodeiv2_control ), intent(in   ) :: c                                    
    type   (fodeiv2_step    ), intent(in   ) :: s                                    
    real   (kind=fgsl_double), intent(in   ) :: dydt                   (:), y0(:), & 
         &                                      yerr                   (:)           
    real   (kind=fgsl_double), intent(inout) :: h                      (:)           
    integer(kind=fgsl_int   )                :: fodeiv2_control_hadjust              
    fodeiv2_control_hadjust = gsl_odeiv2_control_hadjust(c%odeiv2_control, s%odeiv2_step, &
         y0, yerr, dydt, h) 
  end function fodeiv2_control_hadjust
  function fodeiv2_control_name (s)
    type     (fodeiv2_step                   ), intent(in   ) :: s                    
    character(kind=fgsl_char, len=fgsl_strmax)                :: fodeiv2_control_name 
    !
    type     (c_ptr                          )                :: name                 
    name = gsl_odeiv2_control_name(s%odeiv2_step)
    fodeiv2_control_name = fgsl_name(name)
  end function fodeiv2_control_name
  function fodeiv2_evolve_alloc(dim)
    integer(kind=fgsl_size_t), intent(in   ) :: dim                  
    type   (fodeiv2_evolve  )                :: fodeiv2_evolve_alloc 
    fodeiv2_evolve_alloc%odeiv2_evolve = &
         gsl_odeiv2_evolve_alloc(dim)
  end function fodeiv2_evolve_alloc
  function fodeiv2_evolve_apply(e, con, step, dydt, t, t1, h, y) 
    type   (fodeiv2_evolve  ), intent(inout) :: e                             
    type   (fodeiv2_control ), intent(inout) :: con                           
    type   (fodeiv2_step    ), intent(inout) :: step                          
    type   (fodeiv2_system  ), intent(in   ) :: dydt                          
    real   (kind=fgsl_double), intent(inout) :: h                   , t, y(:) 
    real   (kind=fgsl_double), intent(in   ) :: t1                            
    integer(kind=fgsl_int   )                :: fodeiv2_evolve_apply          
    !    write(6, *) 'Start evolving'
    fodeiv2_evolve_apply = gsl_odeiv2_evolve_apply(e%odeiv2_evolve, &
         con%odeiv2_control, step%odeiv2_step, dydt%odeiv2_system, &
         t, t1, h, y)
  end function fodeiv2_evolve_apply
  function fodeiv2_evolve_reset(s)
    type   (fodeiv2_evolve), intent(inout) :: s                    
    integer(kind=c_int    )                :: fodeiv2_evolve_reset 
    fodeiv2_evolve_reset = gsl_odeiv2_evolve_reset(s%odeiv2_evolve)
  end function fodeiv2_evolve_reset
  subroutine fodeiv2_evolve_free(s)
    type(fodeiv2_evolve), intent(inout) :: s 
    call gsl_odeiv2_evolve_free(s%odeiv2_evolve)
  end subroutine fodeiv2_evolve_free
  function fodeiv2_evolve_status(s)
    type   (fodeiv2_evolve), intent(in   ) :: s                     
    logical                                :: fodeiv2_evolve_status 
    fodeiv2_evolve_status = .true.
    if (.not. c_associated(s%odeiv2_evolve)) &
         fodeiv2_evolve_status = .false.
  end function fodeiv2_evolve_status
  function fodeiv2_control_status(s)
    type   (fodeiv2_control), intent(in   ) :: s                      
    logical                                 :: fodeiv2_control_status 
    fodeiv2_control_status = .true.
    if (.not. c_associated(s%odeiv2_control)) &
         fodeiv2_control_status = .false.
  end function fodeiv2_control_status
  function fodeiv2_step_status(s)
    type   (fodeiv2_step), intent(in   ) :: s                   
    logical                              :: fodeiv2_step_status 
    fodeiv2_step_status = .true.
    if (.not. c_associated(s%odeiv2_step)) &
         fodeiv2_step_status = .false.
  end function fodeiv2_step_status
  function fodeiv2_system_status(s)
    type   (fodeiv2_system), intent(in   ) :: s                     
    logical                                :: fodeiv2_system_status 
    fodeiv2_system_status = .true.
    if (.not. c_associated(s%odeiv2_system)) &
         fodeiv2_system_status = .false.
  end function fodeiv2_system_status
  function fodeiv2_driver_apply(d, t, t1, y &
#ifdef PROFILE
       &                        ,sa         &
#endif
       &                       )
    type(fodeiv2_driver  ), intent(in   ) :: d     
    real(kind=fgsl_double), intent(inout) :: t     
    real(kind=fgsl_double), intent(in   ) :: t1    
    real(kind=fgsl_double), intent(inout) :: y (:) 
#ifdef PROFILE
    type(c_funptr), intent(in   ) :: sa 
#endif
    integer(kind=fgsl_int) :: fodeiv2_driver_apply 
    
    fodeiv2_driver_apply = gsl_odeiv2_driver_apply(d%odeiv2_driver, t, t1, y &
#ifdef PROFILE
         & ,sa &
#endif
         & )
  end function fodeiv2_driver_apply
  function fodeiv2_driver_reset(d)
    type   (fodeiv2_driver), intent(inout) :: d                    
    integer(kind=fgsl_int )                :: fodeiv2_driver_reset 
    fodeiv2_driver_reset = gsl_odeiv2_driver_reset(d%odeiv2_driver)
  end function fodeiv2_driver_reset
    subroutine fodeiv2_driver_free(d)
    type(fodeiv2_driver), intent(inout) :: d 
    call gsl_odeiv2_driver_free(d%odeiv2_driver)
  end subroutine fodeiv2_driver_free
  function fodeiv2_driver_alloc_standard_new(system, t, hstart, eps_abs, eps_rel, a_y, a_dydt)
    real  (kind=fgsl_double ), intent(in   ) :: a_dydt                           , a_y    , & 
         &                                      eps_abs                          , eps_rel, & 
         &                                      hstart                                        
    type  (fodeiv2_step_type), intent(in   ) :: t                                             
    type  (fodeiv2_system   ), intent(in   ) :: system                                        
    type  (fodeiv2_driver   )                :: fodeiv2_driver_alloc_standard_new             
    type  (c_ptr            )                :: step_type                                     
    step_type = gsl_odeiv2_aux_odeiv_step_alloc(t%which)
    fodeiv2_driver_alloc_standard_new%odeiv2_driver = &
         gsl_odeiv2_driver_alloc_standard_new(system%odeiv2_system, step_type, hstart, eps_abs, eps_rel, a_y, a_dydt)
  end function fodeiv2_driver_alloc_standard_new
  function fodeiv2_driver_alloc_y_new(system, t, hstart, eps_abs, eps_rel)
    real  (kind=fgsl_double ), intent(in   ) :: eps_abs                   , eps_rel, & 
         &                                      hstart                                 
    type  (fodeiv2_step_type), intent(in   ) :: t                                      
    type  (fodeiv2_system   ), intent(in   ) :: system                                 
    type  (fodeiv2_driver   )                :: fodeiv2_driver_alloc_y_new             
    type  (c_ptr            )                :: step_type                              
    step_type = gsl_odeiv2_aux_odeiv_step_alloc(t%which)
    fodeiv2_driver_alloc_y_new%odeiv2_driver = &
         gsl_odeiv2_driver_alloc_y_new(system%odeiv2_system, step_type, hstart, eps_abs, eps_rel)
  end function fodeiv2_driver_alloc_y_new
  function fodeiv2_driver_alloc_yp_new(system, t, hstart, eps_abs, eps_rel)
    real  (kind=fgsl_double ), intent(in   ) :: eps_abs                    , eps_rel, & 
         &                                      hstart                                  
    type  (fodeiv2_step_type), intent(in   ) :: t                                       
    type  (fodeiv2_system   ), intent(in   ) :: system                                  
    type  (fodeiv2_driver   )                :: fodeiv2_driver_alloc_yp_new             
    type  (c_ptr            )                :: step_type                               
    step_type = gsl_odeiv2_aux_odeiv_step_alloc(t%which)
    fodeiv2_driver_alloc_yp_new%odeiv2_driver = &
         gsl_odeiv2_driver_alloc_yp_new(system%odeiv2_system, step_type, hstart, eps_abs, eps_rel)
  end function fodeiv2_driver_alloc_yp_new
  function fodeiv2_driver_alloc_scaled_new(system, t, hstart, eps_abs, eps_rel, a_y, a_dydt, scale_abs)
    real  (kind=fgsl_double ), intent(in   ) :: a_dydt                            , a_y    , & 
         &                                      eps_abs                           , eps_rel, & 
         &                                      hstart                                         
    real  (kind=fgsl_double ), intent(in   ) :: scale_abs                      (:)             
    type  (fodeiv2_step_type), intent(in   ) :: t                                              
    type  (fodeiv2_system   ), intent(in   ) :: system                                         
    type  (fodeiv2_driver   )                :: fodeiv2_driver_alloc_scaled_new                
    type  (c_ptr            )                :: step_type                                      
    step_type = gsl_odeiv2_aux_odeiv_step_alloc(t%which)
    fodeiv2_driver_alloc_scaled_new%odeiv2_driver = &
         gsl_odeiv2_driver_alloc_scaled_new(system%odeiv2_system, step_type, hstart, eps_abs, eps_rel, a_y, a_dydt, scale_abs)
  end function fodeiv2_driver_alloc_scaled_new
  function fodeiv2_driver_set_hmin(d, hmin)
    real   (kind=fgsl_double), intent(in   ) :: hmin                    
    type   (fodeiv2_driver  ), intent(inout) :: d                       
    integer(kind=c_int      )                :: fodeiv2_driver_set_hmin 
     fodeiv2_driver_set_hmin= &
         gsl_odeiv2_driver_set_hmin(d%odeiv2_driver, hmin)
   end function fodeiv2_driver_set_hmin
   function fodeiv2_driver_set_hmax(d, hmax)
    real   (kind=fgsl_double), intent(in   ) :: hmax                    
    type   (fodeiv2_driver  ), intent(inout) :: d                       
    integer(kind=c_int      )                :: fodeiv2_driver_set_hmax 
     fodeiv2_driver_set_hmax= &
         gsl_odeiv2_driver_set_hmax(d%odeiv2_driver, hmax)
   end function fodeiv2_driver_set_hmax
   function fodeiv2_driver_set_nmax(d, nmax)
   integer(kind=c_long   ), intent(in   ) :: nmax                    
   type   (fodeiv2_driver), intent(inout) :: d                       
   integer(kind=c_int    )                :: fodeiv2_driver_set_nmax 
     fodeiv2_driver_set_nmax= &
         gsl_odeiv2_driver_set_nmax(d%odeiv2_driver, nmax)
   end function fodeiv2_driver_set_nmax
  function fodeiv2_driver_status(d)
    type   (fodeiv2_driver), intent(in   ) :: d                     
    logical                                :: fodeiv2_driver_status 
    fodeiv2_driver_status = .true.
    if (.not. c_associated(d%odeiv2_driver)) &
         fodeiv2_driver_status = .false.
  end function fodeiv2_driver_status
end module fodeiv2

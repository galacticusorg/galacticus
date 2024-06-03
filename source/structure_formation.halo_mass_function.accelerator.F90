!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
Implements a dark matter halo mass function class which accelerates another mass function using tabulation.
  !!}

  use :: Numerical_Interpolation, only : interpolator

  !![
  <haloMassFunction name="haloMassFunctionAccelerator">
   <description>
    A dark matter halo mass function class which accelerates another mass function using tabulation.
   </description>
  </haloMassFunction>
  !!]
  type, extends(haloMassFunctionClass) :: haloMassFunctionAccelerator
     !!{
     A dark matter halo mass function class which accelerates another mass function using tabulation.
     !!}
     private
     class           (haloMassFunctionClass), pointer                   :: haloMassFunction_ => null()
     type            (interpolator         ), allocatable               :: massFunction_              , massFunctionIntegrated_, &
          &                                                                massFraction_
     double precision                       , allocatable, dimension(:) :: mass                       , massFunction           , &
          &                                                                slope
     double precision                                                   :: massMinimum                , massMaximum            , &
          &                                                                time
   contains
     !![
     <methods>
       <method method="tabulate" description="Tabualte the mass function."/>
     </methods>
     !!]
     final     ::                 acceleratorDestructor
     procedure :: differential => acceleratorDifferential
     procedure :: integrated   => acceleratorIntegrated
     procedure :: massFraction => acceleratorMassFraction
     procedure :: tabulate     => acceleratorTabulate
  end type haloMassFunctionAccelerator

  interface haloMassFunctionAccelerator
     !!{
     Constructors for the {\normalfont \ttfamily accelerator} halo mass function class.
     !!}
     module procedure acceleratorConstructorParameters
     module procedure acceleratorConstructorInternal
  end interface haloMassFunctionAccelerator

  ! The number of points per octave of halo mass at which to tabulate the halo mass function.
  double precision, parameter :: countPointsPerOctave=10.0d0
  
contains

  function acceleratorConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily accelerator} halo mass function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (haloMassFunctionAccelerator)                :: self
    type (inputParameters            ), intent(inout) :: parameters
    class(cosmologyParametersClass   ), pointer       :: cosmologyParameters_
    class(haloMassFunctionClass      ), pointer       :: haloMassFunction_

    !![
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="haloMassFunction"    name="haloMassFunction_"    source="parameters"/>
    !!]
    self=haloMassFunctionAccelerator(cosmologyParameters_,haloMassFunction_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="haloMassFunction_"   />
    !!]
    return
  end function acceleratorConstructorParameters

  function acceleratorConstructorInternal(cosmologyParameters_,haloMassFunction_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily accelerator} halo mass function class.
    !!}
    implicit none
    type (haloMassFunctionAccelerator)                        :: self
    class(cosmologyParametersClass   ), target, intent(in   ) :: cosmologyParameters_
    class(haloMassFunctionClass      ), target, intent(in   ) :: haloMassFunction_
    !![
    <constructorAssign variables="*cosmologyParameters_, *haloMassFunction_"/>
    !!]

    self%massMinimum=+huge(0.0d0)
    self%massMaximum=-huge(0.0d0)
    self%time       =-huge(0.0d0)
    return
  end function acceleratorConstructorInternal

  subroutine acceleratorDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily accelerator} halo mass function class.
    !!}
    implicit none
    type(haloMassFunctionAccelerator), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    <objectDestructor name="self%haloMassFunction_"   />
    !!]
    return
  end subroutine acceleratorDestructor

  double precision function acceleratorDifferential(self,time,mass,node) result(massFunction)
    !!{
    Return the differential halo mass function at the given time and mass.
    !!}
    implicit none
    class           (haloMassFunctionAccelerator), intent(inout), target   :: self
    double precision                             , intent(in   )           :: time, mass
    type            (treeNode                   ), intent(inout), optional :: node

    call self%tabulate(time,mass,mass,node)
    massFunction=+exp(self%massFunction_%interpolate(log(mass)))
    return
  end function acceleratorDifferential

  double precision function acceleratorIntegrated(self,time,massLow,massHigh,node,status) result(massFunction)
    !!{
    Return the integrated halo mass function at the given time and mass.
    !!}
    use :: Error, only : errorStatusSuccess
    implicit none
    class           (haloMassFunctionAccelerator), intent(inout), target           :: self
    double precision                             , intent(in   )                   :: time    , massLow, &
         &                                                                            massHigh
    type            (treeNode                   ), intent(inout), target, optional :: node 
    integer                                      , intent(  out)        , optional :: status

    call self%tabulate(time,massLow,massHigh,node)
    massFunction=-exp(self%massFunctionIntegrated_%interpolate(log(massHigh))) &
         &       +exp(self%massFunctionIntegrated_%interpolate(log(massLow )))
    if (present(status)) status=errorStatusSuccess
    return
  end function acceleratorIntegrated
  
  double precision function acceleratorMassFraction(self,time,massLow,massHigh,node) result(massFraction)
    !!{
    Return the integrated halo mass fraction at the given time and mass.
    !!}
    implicit none
    class           (haloMassFunctionAccelerator), intent(inout), target           :: self
    double precision                             , intent(in   )                   :: time    , massLow, &
         &                                                                            massHigh
    type            (treeNode                   ), intent(inout), target, optional :: node 

    call self%tabulate(time,massLow,massHigh,node)
    massFraction=-exp(self%massFraction_%interpolate(log(massHigh))) &
         &       +exp(self%massFraction_%interpolate(log(massLow )))
    return
  end function acceleratorMassFraction

  subroutine acceleratorTabulate(self,time,massLow,massHigh,node)
    !!{
    Tabulate the mass function.
    !!}
    use, intrinsic :: ISO_C_Binding          , only : c_size_t
    use            :: Numerical_Interpolation, only : gsl_interp_linear
    use            :: Numerical_Ranges       , only : Make_Range                  , rangeTypeLogarithmic
    use            :: Table_Labels           , only : extrapolationTypeExtrapolate
    implicit none
    class           (haloMassFunctionAccelerator), intent(inout), target                 :: self
    double precision                             , intent(in   )                         :: time                                                       , massLow     , &
         &                                                                                  massHigh
    type            (treeNode                   ), intent(inout), target      , optional :: node 
    double precision                             , allocatable  , dimension(:)           :: mass                                                       , massFunction, &
         &                                                                                  massFunctionIntegrated                                     , massFraction
    double precision                             , parameter                             :: massRatio             =exp(log(2.0d0)/countPointsPerOctave)
    double precision                             , parameter                             :: massFunctionTiny      =1.0d-100
    double precision                             , parameter                             :: massFunctionHuge      =1.0d+100
    double precision                                                                     :: massMinimum                                                , massMaximum , &
         &                                                                                  slope
    integer         (c_size_t                   )                                        :: countMasses                                                , iMinimum    , &
         &                                                                                  iMaximum                                                   , i

    if (time /= self%time) then
       self%massMinimum=+huge(0.0d0)
       self%massMaximum=-huge(0.0d0)
       self%time       =+time
    end if
    if     (                              &
         &   massLow  >= self%massMinimum &
         &  .and.                         &
         &   massHigh <= self%massMaximum &
         & ) return
    ! Set an initial range of radii that brackets the requested radii.
    massMinimum=0.5d0*massLow
    massMaximum=2.0d0*massHigh
    ! Round to the nearest factor of 2.
    massMinimum=2.0d0**(floor  (log(massMinimum)/log(2.0d0))+0)
    massMaximum=2.0d0**(ceiling(log(massMaximum)/log(2.0d0))+1)
    !! Expand to encompass any pre-existing range.
    if (allocated(self%mass)) then
       massMinimum=min(massMinimum,self%massMinimum)
       massMaximum=max(massMaximum,self%massMaximum)
    end if
    !! Construct arrays.
    countMasses=nint(log(massMaximum/massMinimum)/log(2.0d0)*countPointsPerOctave+1.0d0)
    allocate(mass                  (countMasses))
    allocate(massFunction          (countMasses))
    allocate(massFraction          (countMasses))
    allocate(massFunctionIntegrated(countMasses))
    mass=Make_Range(massMinimum,massMaximum,int(countMasses),rangeTypeLogarithmic)
    ! Copy in any usable results from any previous solution.
    !! Assume by default that no previous solutions are usable.
    iMinimum=+huge(0_c_size_t)
    iMaximum=-huge(0_c_size_t)
    !! Check that a pre-existing solution exists.
    if (allocated(self%mass)) then
       iMinimum=nint(log(self%massMinimum/massMinimum)/log(2.0d0)*countPointsPerOctave)+1_c_size_t
       iMaximum=nint(log(self%massMaximum/massMinimum)/log(2.0d0)*countPointsPerOctave)+1_c_size_t
       massFunction(iMinimum:iMaximum)=self%massFunction
    end if
    ! Compute the mass function where old results were unavailable.
    do i=1,countMasses
       ! Skip cases for which we have a pre-existing solution.
       if (i >= iMinimum .and. i <= iMaximum) cycle
       ! Evaluate the mass function.
       massFunction(i)=min(max(self%haloMassFunction_%differential(time,mass(i),node),massFunctionTiny),massFunctionHuge)
    end do
    ! Compute the integrated mass function.
    massFraction          (countMasses)=massFunctionTiny
    massFunctionIntegrated(countMasses)=massFunctionTiny
    do i=countMasses-1,1,-1
       slope  =+log(                   &
            &       +massFunction(i+1) &
            &       /massFunction(i  ) &
            &      )                   &
            &  /log(                   &
            &       +massRatio         &
            &      )
       if (slope == -1.0d0) then
          massFunctionIntegrated(i)=massFunctionIntegrated(i+1)+massFunction(i)*mass(i)   *log(massRatio                     )
       else
          massFunctionIntegrated(i)=massFunctionIntegrated(i+1)+massFunction(i)*mass(i)   *   (massRatio**(slope+1.0d0)-1.0d0) &
               &                                                                          /               (slope+1.0d0)
       end if
       if (slope == -2.0d0) then
          massFraction          (i)=massFraction          (i+1)+massFunction(i)*mass(i)**2*log(massRatio                     )
       else
          massFraction          (i)=massFraction          (i+1)+massFunction(i)*mass(i)**2*   (massRatio**(slope+2.0d0)-1.0d0) &
               &                                                                          /               (slope+2.0d0)
       end if
    end do
    massFraction=+                          massFraction      &
         &       /self%cosmologyParameters_%densityCritical() &
         &       /self%cosmologyParameters_%OmegaMatter    ()
    ! Build the interpolator.
    if (allocated(self%massFunction_          )) deallocate(self%massFunction_          )
    if (allocated(self%massFraction_          )) deallocate(self%massFraction_          )
    if (allocated(self%massFunctionIntegrated_)) deallocate(self%massFunctionIntegrated_)
    allocate(self%massFunction_          )
    allocate(self%massFraction_          )
    allocate(self%massFunctionIntegrated_)
    self%massFunction_          =interpolator(log(mass),log(massFunction          ),interpolationType=gsl_interp_linear,extrapolationType=extrapolationTypeExtrapolate)
    self%massFraction_          =interpolator(log(mass),log(massFraction          ),interpolationType=gsl_interp_linear,extrapolationType=extrapolationTypeExtrapolate)
    self%massFunctionIntegrated_=interpolator(log(mass),log(massFunctionIntegrated),interpolationType=gsl_interp_linear,extrapolationType=extrapolationTypeExtrapolate)
    ! Store the current results for future re-use.
    if (allocated(self%mass        )) deallocate(self%mass        )
    if (allocated(self%massFunction)) deallocate(self%massFunction)
    allocate(self%mass        (countMasses))
    allocate(self%massFunction(countMasses))
    self%mass        =mass
    self%massFunction=massFunction
    self%massMinimum =massMinimum
    self%massMaximum =massMaximum
    return
  end subroutine acceleratorTabulate


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
  Implementation of a exponential variogram model.
  !!}

  !![
  <variogram name="variogramExponential">
   <description>
   An exponential variogram model in which the variogram is given by:
   \begin{equation}
   C(r) = C_0 + C_1 \left[ 1 - \exp\left( - \frac{r}{C_\mathrm{r}} \right) \right],
   \end{equation}
   where $r$ is the separation, and $(C_0,C_1,C_\mathrm{r})$ are parameters of the model.
   </description>
  </variogram>
  !!]
  type, extends(variogramClass) :: variogramExponential
     !!{
     Implementation of a exponential variogram model.
     !!}
     private
     logical                                             :: assumeZeroVarianceAtZeroLag
     integer         (c_size_t                         ) :: countParameters_
     double precision                                    :: C0                         , C1, &
          &                                                 CR
     type            (enumerationVariogramFitOptionType) :: variogramFitOption
   contains
     procedure :: fit               => exponentialFit
     procedure :: countParameters   => exponentialCountParameters
     procedure :: modelInitialGuess => exponentialModelInitialGuess
     procedure :: modelF            => exponentialModelF
     procedure :: modeldF           => exponentialModeldF
     procedure :: modelFdF          => exponentialModelFdF
     procedure :: variogram         => exponentialVariogram
     procedure :: correlation       => exponentialCorrelation
  end type variogramExponential

  interface variogramExponential
     !!{
     Constructors for the \refClass{variogramExponential} posterior sampling likelihood class.
     !!}
     module procedure exponentialConstructorParameters
     module procedure exponentialConstructorInternal
  end interface variogramExponential

contains

  function exponentialConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{variogramExponential} variogram class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (variogramExponential)                :: self
    type   (inputParameters     ), intent(inout) :: parameters
    logical                                      :: assumeZeroVarianceAtZeroLag
    type   (varying_string      )                :: variogramFitOption

    !![
    <inputParameter>
      <name>variogramFitOption</name>
      <description>Option controlling how the variogram is to be fit.</description>
      <defaultValue>var_str('median')</defaultValue>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>assumeZeroVarianceAtZeroLag</name>
      <description>      
       If true, the variogram model is forced to go to zero for states with zero separation (as expected if the likelihood model
       being emulated is fully deterministic). Otherwise, the variance at zero separation is treated as a free parameter.
      </description>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=variogramExponential(enumerationVariogramFitOptionEncode(char(variogramFitOption),includesPrefix=.false.),assumeZeroVarianceAtZeroLag)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function exponentialConstructorParameters

  function exponentialConstructorInternal(variogramFitOption,assumeZeroVarianceAtZeroLag) result(self)
    !!{
    Constructor for the \refClass{variogramExponential} variogram class.
    !!}
    implicit none
    type   (variogramExponential             )                :: self
    type   (enumerationVariogramFitOptionType), intent(in   ) :: variogramFitOption
    logical                                   , intent(in   ) :: assumeZeroVarianceAtZeroLag
    !![
    <constructorAssign variables="variogramFitOption, assumeZeroVarianceAtZeroLag"/>
    !!]
    
    if (self%assumeZeroVarianceAtZeroLag) then
       self%countParameters_=2_c_size_t
    else
       self%countParameters_=3_c_size_t
    end if
    return
  end function exponentialConstructorInternal

  function exponentialCountParameters(self)
    !!{
    Return the number of parameters in the exponential variogram model.
    !!}
    implicit none
    integer(c_size_t            )                :: exponentialCountParameters
    class  (variogramExponential), intent(inout) :: self
    
    exponentialCountParameters=self%countParameters_
    return
  end function exponentialCountParameters

  double precision function exponentialVariogram(self,separation)
    !!{
    Compute the variogram at the given {\normalfont \ttfamily separation}. Here we use the ``exponential model'':
    \begin{equation}
      \gamma(h) =C_0 + C_1 [ 1 - \exp\left(-\frac{h}{a}\right)]
    \end{equation}
    where $h$ is the separation, and $C_0$, $C_1$, and $R_\mathrm{C}$ are parameters of the model.
    !!}
    implicit none
    class           (variogramExponential), intent(inout)           :: self
    double precision                      , intent(in   ), optional :: separation

    if (present(separation)) then
       exponentialVariogram=+self%C0                    &
            &               +self%C1                    &
            &               *(                          &
            &                 +1.0d0                    &
            &                 -exp(-separation/self%CR) &
            &                )
    else
       exponentialVariogram=+self%C0                    &
            &               +self%C1
    end if
    return
  end function exponentialVariogram

  double precision function exponentialCorrelation(self,separation)
    !!{
    Compute correlation of the variogram model. This is just the variance at maximum separation minus the variogram.
    !!}
    implicit none
    class           (variogramExponential), intent(inout) :: self
    double precision                      , intent(in   ) :: separation

    exponentialCorrelation=+self%C0                    &
         &                 +self%C1                    &
         &                 -self%variogram(separation)
    return
  end function exponentialCorrelation

  function exponentialModelInitialGuess(self,separations,semiVariances) result(C)
    !!{
    Provide an initial guess for parameters of the exponential variogram model.
    !!}
    implicit none
    double precision                      , allocatable  , dimension(:) :: C
    class           (variogramExponential), intent(inout)               :: self
    double precision                      , intent(in   ), dimension(:) :: separations, semiVariances

    allocate(C(self%countParameters_))
    if (self%assumeZeroVarianceAtZeroLag) then
       C=[                 semiVariances(size(semiVariances)),separations(size(separations)/2)]
    else
       C=[semiVariances(1),semiVariances(size(semiVariances)),separations(size(separations)/2)]
    end if
    ! Convert to logarithmic form.
    C=log(C)
    return
  end function exponentialModelInitialGuess

  subroutine exponentialFit(self,separations,semiVariances)
    !!{
    Compute best fit coefficients for the exponential variogram model.
    !!}
    implicit none
    class           (variogramExponential), intent(inout)                            :: self
    double precision                      , intent(in   )             , dimension(:) :: separations, semiVariances
    double precision                                     , allocatable, dimension(:) :: C
    
    call self%fitGeneric(self%variogramFitOption,separations,semiVariances,C)
    if (self%assumeZeroVarianceAtZeroLag) then
       self%C0=0.0d0
       self%C1=exp(C(1))
       self%CR=exp(C(2))
    else
       self%C0=exp(C(1))
       self%C1=exp(C(2))
       self%CR=exp(C(3))       
    end if
    self%C0=self%C0*self%semiVarianceNormalization
    self%C1=self%C1*self%semiVarianceNormalization
    self%CR=self%CR*self%  separationNormalization
    return
  end subroutine exponentialFit
  
  double precision function exponentialModelF(self,C,separations,semiVariances)
    !!{
    Function to be minimized when fitting the variogram.
    !!}
    implicit none
    class           (variogramExponential), intent(inout)                     :: self
    double precision                      , intent(in   ), dimension(     : ) :: C            , separations, &
         &                                                                       semiVariances
    double precision                                     , dimension(size(C)) :: C_
    double precision                                                          :: C0           , C1         , &
         &                                                                       CR

    C_=exp(C)
    if (self%assumeZeroVarianceAtZeroLag) then
       C0=0.0d0
       C1=C_(1)
       CR=C_(2)
    else
       C0=C_(1)
       C1=C_(2)
       CR=C_(3)
    end if
    exponentialModelF=sum(                               &
         &                +(                             &
         &                  +log(                        &
         &                       +C0                     &
         &                       +C1                     &
         &                       *(                      &
         &                         +1.0d0                &
         &                         -exp(-separations/CR) &
         &                        )                      &
         &                      )                        &
         &                  -log(                        &
         &                       +semiVariances          &
         &                      )                        &
         &                 )**2                          &
         &               )
    return
  end function exponentialModelF

  function exponentialModeldF(self,C,separations,semiVariances) result(dfdC)
    !!{
    Derivatives of the function to be minimized when fitting the variogram.
    !!}
    implicit none
    class           (variogramExponential), intent(inout)                               :: self
    double precision                      , intent(in   ), dimension(               : ) :: C            , separations, &
         &                                                                                 semivariances
    double precision                      , allocatable  , dimension(               : ) :: dfdC
    double precision                                     , dimension(size(          C)) :: C_
    double precision                                     , dimension(size(separations)) :: f
    double precision                                                                    :: C0           , C1         , &
         &                                                                                 CR
    integer                                                                             :: i1           , iR

    allocate(dfdC(size(C)))
    C_=exp(C)
    if (self%assumeZeroVarianceAtZeroLag) then
       C0=0.0d0
       i1=1
       iR=2
    else
       C0=C_(1)
       i1=2
       iR=3
    end if
    C1=C_(i1)
    CR=C_(iR)
    f      =+(                             &
         &    +log(                        &
         &         +C0                     &
         &         +C1                     &
         &         *(                      &
         &           +1.0d0                &
         &           -exp(-separations/CR) &
         &          )                      &
         &        )                        &
         &    -log(                        &
         &         +semiVariances          &
         &        )                        &
         &   )                             &
         &  /(                             &
         &         +C0                     &
         &         +C1                     &
         &         *(                      &
         &           +1.0d0                &
         &           -exp(-separations/CR) &
         &          )                      &
         &   )
    if (.not.self%assumeZeroVarianceAtZeroLag) &
         & dfdC(1)=+2.0d0*sum(f                                                   )
    dfdC(i1)      =+2.0d0*sum(f   *(+1.0d0-exp(-separations/CR)                  ))
    dfdC(iR)      =+2.0d0*sum(f*C1*(      -exp(-separations/CR)*separations/CR**2))
    dfdC          =+dfdC                                                            &
         &         *C_
    where(abs(dfdC) < 1.0d-30)
       dfdC=1.0d-30
    end where
    return
  end function exponentialModeldF

  subroutine exponentialModelFdF(self,C,separations,semiVariances,f,dfdC)
    !!{
    Computes both function and derivatives to be minimized when fitting the variogram.
    !!}
    implicit none
    class           (variogramExponential), intent(inout)                            :: self
    double precision                      , intent(in   )             , dimension(:) :: C            , separations, &
         &                                                                              semivariances
    double precision                      , intent(  out)                            :: f
    double precision                      , intent(  out), allocatable, dimension(:) :: dfdC

    f   =self%modelF (C,separations,semiVariances)
    dfdC=self%modeldF(C,separations,semiVariances)
    return
  end subroutine exponentialModelFdF

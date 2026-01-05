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
  Implements an excursion set barrier class which remaps another class using the \cite{sheth_ellipsoidal_2001} ellipsoidal collapse parameterization.
  !!}

  !![
  <excursionSetBarrier name="excursionSetBarrierRemapShethMoTormen">
   <description>
    An excursion set barrier class which remaps another class using the \cite{sheth_ellipsoidal_2001} ellipsoidal collapse
    parameterization:
    \begin{equation}
     B(S) \rightarrow \sqrt{A} B(S) \left(1 + b \left[ {S \over A B^2(S)}\right]^c\right),
    \end{equation}
    where $A=0.707$, $b=0.5$, and $c=0.6$.
   </description>
  </excursionSetBarrier>
  !!]
  type, extends(excursionSetBarrierClass) :: excursionSetBarrierRemapShethMoTormen
     !!{
     An excursion set barrier class which remaps another class using the \cite{sheth_ellipsoidal_2001} ellipsoidal collapse parameterization.
     !!}
     private
     class           (excursionSetBarrierClass        ), pointer :: excursionSetBarrier_ => null()
     double precision                                            :: a                             , b, &
          &                                                         c
     type            (enumerationExcursionSetRemapType)          :: applyTo
   contains
     final     ::                    remapShethMoTormenDestructor
     procedure :: barrier         => remapShethMoTormenBarrier
     procedure :: barrierGradient => remapShethMoTormenBarrierGradient
  end type excursionSetBarrierRemapShethMoTormen

  interface excursionSetBarrierRemapShethMoTormen
     !!{
     Constructors for the remap scale excursion set barrier class.
     !!}
     module procedure remapShethMoTormenConstructorParameters
     module procedure remapShethMoTormenConstructorInternal
  end interface excursionSetBarrierRemapShethMoTormen

contains

  function remapShethMoTormenConstructorParameters(parameters) result(self)
    !!{
    Constructor for the critical overdensity excursion set class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (excursionSetBarrierRemapShethMoTormen)                :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    class           (excursionSetBarrierClass             ), pointer       :: excursionSetBarrier_
    double precision                                                       :: a                   , b, &
         &                                                                    c
    type            (varying_string                       )                :: applyTo
    
    ! Check and read parameters.
    !![
    <inputParameter>
      <name>a</name>
      <source>parameters</source>
      <variable>a</variable>
      <defaultValue>0.707d0</defaultValue>
      <description>The parameter $a$ in the \cite{sheth_ellipsoidal_2001} ellipsoidal collapse excursion set barrier remapping.</description>
    </inputParameter>
    <inputParameter>
      <name>b</name>
      <source>parameters</source>
      <variable>b</variable>
      <defaultValue>0.500d0</defaultValue>
      <description>The parameter $b$ in the \cite{sheth_ellipsoidal_2001} ellipsoidal collapse excursion set barrier remapping.</description>
    </inputParameter>
    <inputParameter>
      <name>c</name>
      <source>parameters</source>
      <variable>c</variable>
      <defaultValue>0.600d0</defaultValue>
      <description>The parameter $c$ in the \cite{sheth_ellipsoidal_2001} ellipsoidal collapse excursion set barrier remapping.</description>
    </inputParameter>
    <inputParameter>
      <name>applyTo</name>
      <source>parameters</source>
      <defaultValue>var_str('nonRates')</defaultValue>
      <description>Specifies whether rescaling is to be applied to the barrier when used for rate calculation, for other calculations, or both.</description>
    </inputParameter>
    <objectBuilder class="excursionSetBarrier" name="excursionSetBarrier_" source="parameters"/>
    !!]
    self=excursionSetBarrierRemapShethMoTormen(a,b,c,enumerationExcursionSetRemapEncode(applyTo,includesPrefix=.false.),excursionSetBarrier_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="excursionSetBarrier_"/>
    !!]
    return
  end function remapShethMoTormenConstructorParameters

  function remapShethMoTormenConstructorInternal(a,b,c,applyTo,excursionSetBarrier_) result(self)
    !!{
    Internal constructor for the critical overdensity excursion set class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (excursionSetBarrierRemapShethMoTormen)                        :: self
    class           (excursionSetBarrierClass             ), intent(in   ), target :: excursionSetBarrier_
    double precision                                       , intent(in   )         :: a                   , b, &
         &                                                                            c
    type            (enumerationExcursionSetRemapType     ), intent(in   )         :: applyTo
    !![
    <constructorAssign variables="a, b, c, applyTo, *excursionSetBarrier_"/>
    !!]

    if (.not.enumerationExcursionSetRemapIsValid(applyTo)) call Error_Report('applyTo is invalid'//{introspection:location})
    return
  end function remapShethMoTormenConstructorInternal

  subroutine remapShethMoTormenDestructor(self)
    !!{
    Destructor for the critical overdensity excursion set barrier class.
    !!}
    implicit none
    type(excursionSetBarrierRemapShethMoTormen), intent(inout) :: self

    !![
    <objectDestructor name="self%excursionSetBarrier_"/>
    !!]
    return
  end subroutine remapShethMoTormenDestructor

  double precision function remapShethMoTormenBarrier(self,variance,time,node,rateCompute)
    !!{
    Return the excursion set barrier at the given variance and time.
    !!}
    implicit none
    class           (excursionSetBarrierRemapShethMoTormen), intent(inout) :: self
    double precision                                       , intent(in   ) :: variance   , time
    type            (treeNode                             ), intent(inout) :: node
    logical                                                , intent(in   ) :: rateCompute

    remapShethMoTormenBarrier=self%excursionSetBarrier_%barrier(variance,time,node,rateCompute)
    if     (                                                                    &
         &    self%applyTo == excursionSetRemapBoth                             &
         &  .or.                                                                &
         &   (self%applyTo == excursionSetRemapRates    .and.      rateCompute) &
         &  .or.                                                                &
         &   (self%applyTo == excursionSetRemapNonRates .and. .not.rateCompute) &
         & )                                                                    &
         & remapShethMoTormenBarrier=+sqrt(self%a)                     &
         &                           *remapShethMoTormenBarrier        &
         &                           *(                                &
         &                             +1.0d0                          &
         &                             +self%b                         &
         &                             *(                              &
         &                               +variance                     &
         &                               /self%a                       &
         &                               /remapShethMoTormenBarrier**2 &
         &                              )**self%c                      &
         &                            )
    return
  end function remapShethMoTormenBarrier

  double precision function remapShethMoTormenBarrierGradient(self,variance,time,node,rateCompute)
    !!{
    Return the gradient with respect to variance of the excursion set barrier at the given variance and time.
    !!}
    implicit none
    class           (excursionSetBarrierRemapShethMoTormen), intent(inout) :: self
    double precision                                       , intent(in   ) :: variance   , time
    type            (treeNode                             ), intent(inout) :: node
    logical                                                , intent(in   ) :: rateCompute
    double precision                                                       :: barrier    , barrierGradient

    barrierGradient=self%excursionSetBarrier_%barrierGradient(variance,time,node,rateCompute)
    if     (                                                                    &
         &    self%applyTo == excursionSetRemapBoth                             &
         &  .or.                                                                &
         &   (self%applyTo == excursionSetRemapRates    .and.      rateCompute) &
         &  .or.                                                                &
         &   (self%applyTo == excursionSetRemapNonRates .and. .not.rateCompute) &
         & )                                                                    &
         & then
       if (variance <= 0.0d0) then
          remapShethMoTormenBarrierGradient=0.0d0
       else
          barrier                          =self%excursionSetBarrier_%barrier(variance,time,node,rateCompute)
          remapShethMoTormenBarrierGradient=+sqrt(self%a)*barrierGradient*(                                      &
               &                                                           +1.0d0                                &
               &                                                           +self%b                               &
               &                                                           *(variance/self%a/barrier**2)**self%c &
               &                                                          )                                      &
               &                            +sqrt(self%a)*barrier        *(                                      &
               &                                                           +self%b                               &
               &                                                           *self%c                               &
               &                                                           *(variance/self%a/barrier**2)**self%c &
               &                                                           *(                                    &
               &                                                             +1.0d0                              &
               &                                                             -2.0d0                              &
               &                                                             *variance                           &
               &                                                             *barrierGradient                    &
               &                                                             /barrier                            &
               &                                                            )                                    &
               &                                                           /variance                             &
               &                                                          )
       end if
    else
       remapShethMoTormenBarrierGradient=+barrierGradient
    end if
    return
  end function remapShethMoTormenBarrierGradient

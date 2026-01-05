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

!+ Contributions to this file made by: Daniel McAndrew.

  !!{
  An implementation of the intergalactic medium state class for an internal model of instantaneous and full reionization.
  !!}

  use :: Numerical_Interpolation, only : interpolator

  !![
  <intergalacticMediumState name="intergalacticMediumStateInternal">
   <description>The state of the intergalactic medium is solved for internally.</description>
  </intergalacticMediumState>
  !!]
  type, extends(intergalacticMediumStateClass) :: intergalacticMediumStateInternal
     !!{
     An \gls{igm} state class for an internally consistent model.
     !!}
     double precision              , allocatable, dimension(:  ) :: time            , temperatureIGM  , &
          &                                                         massFiltering   , densityHydrogen1, &
          &                                                         densityHydrogen2, densityHelium1  , &
          &                                                         densityHelium2  , densityHelium3
     type            (interpolator), allocatable                 :: interpolator_
   contains
     final     ::                                internalDestructor
     procedure :: electronFraction            => internalElectronFraction
     procedure :: temperature                 => internalTemperature
     procedure :: neutralHydrogenFraction     => internalNeutralHydrogenFraction
     procedure :: neutralHeliumFraction       => internalNeutralHeliumFraction
     procedure :: singlyIonizedHeliumFraction => internalSinglyIonizedHeliumFraction
     procedure :: autoHook                    => internalAutoHook
  end type intergalacticMediumStateInternal

  interface intergalacticMediumStateInternal
     !!{
     Constructors for the internal intergalactic medium state class.
     !!}
     module procedure internalConstructorParameters
     module procedure internalConstructorInternal
  end interface intergalacticMediumStateInternal

contains

  function internalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{intergalacticMediumStateInternal} \gls{igm} state class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (intergalacticMediumStateInternal)                :: self
    type (inputParameters                 ), intent(inout) :: parameters
    class(cosmologyFunctionsClass         ), pointer       :: cosmologyFunctions_
    class(cosmologyParametersClass        ), pointer       :: cosmologyParameters_

    !![
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !!]
    self=intergalacticMediumStateInternal(cosmologyFunctions_,cosmologyParameters_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="cosmologyParameters_"/>
    !!]
    return
  end function internalConstructorParameters

  function internalConstructorInternal(cosmologyFunctions_,cosmologyParameters_) result(self)
    !!{
    Internal constructor for the \refClass{intergalacticMediumStateInternal} \gls{igm} state class.
    !!}
    implicit none
    type (intergalacticMediumStateInternal)                        :: self
    class(cosmologyFunctionsClass         ), intent(inout), target :: cosmologyFunctions_
    class(cosmologyParametersClass        ), intent(inout), target :: cosmologyParameters_
    !![
    <constructorAssign variables="*cosmologyFunctions_, *cosmologyParameters_"/>
    !!]

    allocate  (self%time            (0))
    allocate  (self%temperatureIGM  (0))
    allocate  (self%densityHydrogen1(0))
    allocate  (self%densityHydrogen2(0))
    allocate  (self%densityHelium1  (0))
    allocate  (self%densityHelium2  (0))
    allocate  (self%densityHelium3  (0))
    allocate  (self%massFiltering   (0))
    deallocate(self%time               )
    deallocate(self%temperatureIGM     )
    deallocate(self%densityHydrogen1   )
    deallocate(self%densityHydrogen2   )
    deallocate(self%densityHelium1     )
    deallocate(self%densityHelium2     )
    deallocate(self%densityHelium3     )
    deallocate(self%massFiltering      )
    return
  end function internalConstructorInternal

  subroutine internalAutoHook(self)
    !!{
    Hook into the internal intergalactic medium state evolver to receive updates.
    !!}
    use :: Events_Hooks, only : intergalacticMediumStateEvolveUpdateEventGlobal
    implicit none
    class(intergalacticMediumStateInternal), intent(inout) :: self

    call intergalacticMediumStateEvolveUpdateEventGlobal%attach(self,internalStateSet,label='intergalacticMediumStateInternal')
    return
  end subroutine internalAutoHook

  subroutine internalDestructor(self)
    !!{
    Destructor for the internal \gls{igm} state class.
    !!}
    use :: Events_Hooks, only : intergalacticMediumStateEvolveUpdateEventGlobal
    implicit none
    type(intergalacticMediumStateInternal), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    <objectDestructor name="self%cosmologyFunctions_" />
    !!]
    if (intergalacticMediumStateEvolveUpdateEventGlobal%isAttached(self,internalStateSet)) call intergalacticMediumStateEvolveUpdateEventGlobal%detach(self,internalStateSet)
    return
  end subroutine internalDestructor

  double precision function internalElectronFraction(self,time)
    !!{
    Return the electron fraction of the \gls{igm} in the internal model.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    class           (intergalacticMediumStateInternal), intent(inout)   :: self
    double precision                                  , intent(in   )   :: time
    double precision                                  , dimension(0:1)  :: h
    integer         (c_size_t                        )                  :: i              , j
    double precision                                                    :: densityHydrogen, densityHelium

    if (size(self%time) > 1) then
       i=self%interpolator_%locate(time)
    else
       densityHydrogen         =self%densityHydrogen1(1)+self%densityHydrogen2(1)
       densityHelium           =self%densityHelium1  (1)+self%densityHelium2  (1)+self%densityHelium3(1)
       internalElectronFraction=+(self%densityHydrogen2(1)+self%densityHelium2(1)+2.0d0*self%densityHelium3(1)) &
            &                   /(densityHydrogen+2.0d0*densityHelium)
       return
    end if
    if (self%time(i+1)-self%time(i) > 0.0d0) then
       call self%interpolator_%linearWeights(time,i,h)
    else
       h(0)=0.0d0
       h(1)=1.0d0
    end if
    internalElectronFraction = 0.0d0
    do j=0,1
       densityHydrogen=self%densityHydrogen1(i+j)+self%densityHydrogen2(i+j)
       densityHelium  =self%densityHelium1  (i+j)+self%densityHelium2  (i+j)+self%densityHelium3(i+j)
       if (densityHydrogen > 0.0d0) internalElectronFraction=+internalElectronFraction           &
            &                                                +h                            (  j) &
            &                                                *(                                  &
            &                                                  +      self%densityHydrogen2(i+j) &
            &                                                  +      self%densityHelium2  (i+j) &
            &                                                  +2.0d0*self%densityHelium3  (i+j) &
            &                                                )                                   &
            &                                                /densityHydrogen
    end do
    internalElectronFraction=max(0.0d0,internalElectronFraction)
    return
  end function internalElectronFraction

  double precision function internalNeutralHydrogenFraction(self,time)
    !!{
    Return the neutral hydrogen fraction of the \gls{igm} in the internal model.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    class           (intergalacticMediumStateInternal), intent(inout)   :: self
    double precision                                  , intent(in   )   :: time
    double precision                                  , dimension(0:1)  :: h
    integer         (c_size_t                        )                  :: i              , j
    double precision                                                    :: densityHydrogen

    if (size(self%time) > 1) then
       i=self%interpolator_%locate(time)
    else
       densityHydrogen                =self%densityHydrogen1(1)+self%densityHydrogen2(1)
       internalNeutralHydrogenFraction=self%densityHydrogen1(1)/densityHydrogen
       return
    end if
    if (self%time(i+1)-self%time(i) > 0.0d0) then
       call self%interpolator_%linearWeights(time,i,h)
    else
       h(0)=0.0d0
       h(1)=1.0d0
    end if
    internalNeutralHydrogenFraction=0.0d0
    do j=0,1
       densityHydrogen=self%densityHydrogen1(i+j)+self%densityHydrogen2(i+j)
       if (densityHydrogen > 0.0d0) internalNeutralHydrogenFraction=+internalNeutralHydrogenFraction &
            &                                                       +h                    (  j)      &
            &                                                       *self%densityHydrogen1(i+j)      &
            &                                                       /densityHydrogen
    end do
    internalNeutralHydrogenFraction=max(0.0d0,min(1.0d0,internalNeutralHydrogenFraction))
    return
  end function internalNeutralHydrogenFraction

  double precision function internalNeutralHeliumFraction(self,time)
    !!{
    Return the neutral helium fraction of the \gls{igm} in the internal model.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    class           (intergalacticMediumStateInternal), intent(inout)  :: self
    double precision                                  , intent(in   )  :: time
    double precision                                  , dimension(0:1) :: h
    integer         (c_size_t                        )                  :: i            , j
    double precision                                                    :: densityHelium

    if (size(self%time) > 1) then
       i=self%interpolator_%locate(time)
    else
       densityHelium                =self%densityHelium1(1)+self%densityHelium2(1)+self%densityHelium3(1)
       internalNeutralHeliumFraction=self%densityHelium1(1)/densityHelium
       return
    end if
    if (self%time(i+1)-self%time(i) > 0.0d0) then
       call self%interpolator_%linearWeights(time,i,h)
    else
       h(0)=0.0d0
       h(1)=1.0d0
    end if
    internalNeutralHeliumFraction=0.0d0
    do j=0,1
       densityHelium=self%densityHelium1(i+j)+self%densityHelium2(i+j)+self%densityHelium3(i+j)
       if (densityHelium > 0.0d0) internalNeutralHeliumFraction=+internalNeutralHeliumFraction &
            &                                                   +h                  (  j)      &
            &                                                   *self%densityHelium1(i+j)      &
            &                                                   /densityHelium
    end do
    internalNeutralHeliumFraction=max(0.0d0,min(1.0d0,internalNeutralHeliumFraction))
    return
  end function internalNeutralHeliumFraction

  double precision function internalSinglyIonizedHeliumFraction(self,time)
    !!{
    Return the singly ionized helium fraction of the \gls{igm} in the internal model.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    class           (intergalacticMediumStateInternal), intent(inout)   :: self
    double precision                                  , intent(in   )   :: time
    double precision                                  , dimension(0:1)  :: h
    integer         (c_size_t                        )                  :: i            , j
    double precision                                                    :: densityHelium

    if (size(self%time) > 1) then
       i=self%interpolator_%locate(time)
    else
       densityHelium                      =self%densityHelium1(1)+self%densityHelium2(1)+self%densityHelium3(1)
       internalSinglyIonizedHeliumFraction=self%densityHelium2(1)/densityHelium
       return
    end if
    if (self%time(i+1)-self%time(i) > 0.0d0) then
       call self%interpolator_%linearWeights(time,i,h)
    else
       h(0)=0.0d0
       h(1)=1.0d0
    end if
    internalSinglyIonizedHeliumFraction=0.0d0
    do j=0,1
       densityHelium=self%densityHelium1(i+j)+self%densityHelium2(i+j)+self%densityHelium3(i+j)
       if (densityHelium > 0.0d0) internalSinglyIonizedHeliumFraction=+internalSinglyIonizedHeliumFraction &
            &                                                         +h                  (  j)            &
            &                                                         *self%densityHelium2(i+j)            &
            &                                                         /densityHelium
    end do
    internalSinglyIonizedHeliumFraction=max(0.0d0,min(1.0d0,internalSinglyIonizedHeliumFraction))
    return
  end function internalSinglyIonizedHeliumFraction

  double precision function internalTemperature(self,time)
    !!{
    Return the temperature of the \gls{igm} in the internal model.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    class           (intergalacticMediumStateInternal), intent(inout)   :: self
    double precision                                  , intent(in   )   :: time
    double precision                                  , dimension(0:1)  :: h
    integer         (c_size_t                        )                  :: i   , j

    if (size(self%time) > 1) then
       i=self%interpolator_%locate(time)
    else
       internalTemperature=self%temperatureIGM(1)
       return
    end if
    if (self%time(i+1)-self%time(i) > 0.0d0) then
       call self%interpolator_%linearWeights(time,i,h)
    else
       h(0)=0.0d0
       h(1)=1.0d0
    end if
    internalTemperature=0.0d0
    do j=0,1
       internalTemperature=internalTemperature+h(j)*self%temperatureIGM(i+j)
    end do
    return
  end function internalTemperature

  subroutine internalStateSet(self,time,densityHydrogen1,densityHydrogen2,densityHelium1,densityHelium2,densityHelium3,temperature,massFiltering)
    !!{
    Set state in the internal intergalactic medium state class.
    !!}
    use :: Error            , only : Error_Report
    implicit none
    class           (*), intent(inout)               :: self
    double precision   , intent(in   ), dimension(:) :: time            , densityHydrogen1, &
         &                                              densityHydrogen2, densityHelium1  , &
         &                                              densityHelium2  , densityHelium3  , &
         &                                              temperature     , massFiltering
    select type (self)
    class is (intergalacticMediumStateInternal)
       if (allocated(self%time            )) deallocate(self%time            )
       if (allocated(self%temperatureIGM  )) deallocate(self%temperatureIGM  )
       if (allocated(self%densityHydrogen1)) deallocate(self%densityHydrogen1)
       if (allocated(self%densityHydrogen2)) deallocate(self%densityHydrogen2)
       if (allocated(self%densityHelium1  )) deallocate(self%densityHelium1  )
       if (allocated(self%densityHelium2  )) deallocate(self%densityHelium2  )
       if (allocated(self%densityHelium3  )) deallocate(self%densityHelium3  )
       if (allocated(self%massFiltering   )) deallocate(self%massFiltering   )
       allocate(self%time            ,mold=time            )
       allocate(self%temperatureIGM  ,mold=temperature     )
       allocate(self%densityHydrogen1,mold=densityHydrogen1)
       allocate(self%densityHydrogen2,mold=densityHydrogen2)
       allocate(self%densityHelium1  ,mold=densityHelium1  )
       allocate(self%densityHelium2  ,mold=densityHelium2  )
       allocate(self%densityHelium3  ,mold=densityHelium3  )
       allocate(self%massFiltering   ,mold=massFiltering   )
       self%time            =time
       self%temperatureIGM  =temperature
       self%densityHydrogen1=densityHydrogen1
       self%densityHydrogen2=densityHydrogen2
       self%densityHelium1  =densityHelium1
       self%densityHelium2  =densityHelium2
       self%densityHelium3  =densityHelium3
       self%massFiltering   =massFiltering
       ! Build interpolator.
       if (allocated(self%interpolator_)) deallocate(self%interpolator_)
       allocate(self%interpolator_)
       self%interpolator_=interpolator(self%time)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine internalStateSet

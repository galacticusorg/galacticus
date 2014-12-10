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

!+ Contributions to this file made by: Daniel McAndrew.

  !% An implementation of the intergalactic medium state class for an internal model of instantaneous and full reionization.

  !# <intergalacticMediumState name="intergalacticMediumStateInternal" defaultThreadPrivate="no">
  !#  <description>The intergalactic medium is assumed to be instantaneously and fully reionized at a fixed redshift, and heated to a fixed temperature.</description>
  !# </intergalacticMediumState>

  type, extends(intergalacticMediumStateClass) :: intergalacticMediumStateInternal
     !% An \gls{igm} state class for an internally consistent model.
     double precision, allocatable, dimension(:  ) :: time            , temperatureIGM  , &
          &                                           massFiltering   , densityHydrogen1, &
          &                                           densityHydrogen2, densityHelium1  , &
          &                                           densityHelium2  , densityHelium3
   contains
     !@ <objectMethods>
     !@   <object>intergalacticMediumStateInternal</object>
     !@   <objectMethod>
     !@     <method>densityH1Set</method>
     !@     <type>\void</type>
     !@     <arguments>\doubleone\ densityHydrogen1\argin</arguments>
     !@     <description>Set the density of neutral hydrogen time series.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>densityH2Set</method>
     !@     <type>\void</type>
     !@     <arguments>\doubleone\ densityHydrogen2\argin</arguments>
     !@     <description>Set the density of ionized hydrogen time series.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>densityHe1Set</method>
     !@     <type>\void</type>
     !@     <arguments>\doubleone\ densityHelium1\argin</arguments>
     !@     <description>Set the density of neutral helium time series.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>densityHe2Set</method>
     !@     <type>\void</type>
     !@     <arguments>\doubleone\ densityHelium2\argin</arguments>
     !@     <description>Set the density of singly-ionized helium time series.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>densityHe3Set</method>
     !@     <type>\void</type>
     !@     <arguments>\doubleone\ densityHelium3\argin</arguments>
     !@     <description>Set the density of doubly-ionized helium time series.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>timeSet</method>
     !@     <type>\void</type>
     !@     <arguments>\doubleone\ times\argin</arguments>
     !@     <description>Set the times to use for all time series.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>temperatureSet</method>
     !@     <type>\void</type>
     !@     <arguments>\doubleone\ temperature\argin</arguments>
     !@     <description>Set the temperature time series.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>massFilteringSet</method>
     !@     <type>\void</type>
     !@     <arguments>\doubleone\ massFiltering\argin</arguments>
     !@     <description>Set the filtering mass time series.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>filteringMAss</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ time\argin</arguments>
     !@     <description>Return the filtering mass at the given {\tt time}.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: electronFraction            => internalElectronFraction
     procedure :: temperature                 => internalTemperature
     procedure :: neutralHydrogenFraction     => internalNeutralHydrogenFraction
     procedure :: neutralHeliumFraction       => internalNeutralHeliumFraction
     procedure :: singlyIonizedHeliumFraction => internalSinglyIonizedHeliumFraction
     procedure :: filteringMass               => internalFilteringMass
     procedure :: timeSet                     => internalTimeSet
     procedure :: temperatureSet              => internalTemperatureSet
     procedure :: densityH1Set                => internalDensityHydrogen1Set
     procedure :: densityH2Set                => internalDensityHydrogen2Set
     procedure :: densityHe1Set               => internalDensityHelium1Set
     procedure :: densityHe2Set               => internalDensityHelium2Set
     procedure :: densityHe3Set               => internalDensityHelium3Set
     procedure :: massFilteringSet            => internalMassFilteringSet
  end type intergalacticMediumStateInternal

  interface intergalacticMediumStateInternal
     !% Constructors for the internal intergalactic medium state class.
     module procedure internalDefaultConstructor
  end interface intergalacticMediumStateInternal

contains

  function internalDefaultConstructor()
    !% Default constructor for the {\tt internal} \gls{igm} state class.
    use Input_Parameters
    implicit none
    type    (intergalacticMediumStateInternal), target  :: internalDefaultConstructor

    allocate  (internalDefaultConstructor%time            (0))
    allocate  (internalDefaultConstructor%temperatureIGM  (0))
    allocate  (internalDefaultConstructor%densityHydrogen1(0))
    allocate  (internalDefaultConstructor%densityHydrogen2(0))
    allocate  (internalDefaultConstructor%densityHelium1  (0))
    allocate  (internalDefaultConstructor%densityHelium2  (0))
    allocate  (internalDefaultConstructor%densityHelium3  (0))
    allocate  (internalDefaultConstructor%massFiltering   (0))
    deallocate(internalDefaultConstructor%time               )
    deallocate(internalDefaultConstructor%temperatureIGM     )
    deallocate(internalDefaultConstructor%densityHydrogen1   )
    deallocate(internalDefaultConstructor%densityHydrogen2   )
    deallocate(internalDefaultConstructor%densityHelium1     )
    deallocate(internalDefaultConstructor%densityHelium2     )
    deallocate(internalDefaultConstructor%densityHelium3     )
    deallocate(internalDefaultConstructor%massFiltering      )
    return
  end function internalDefaultConstructor

  double precision function internalFilteringMass(self,time)
    !% Return the filtering mass of the \gls{igm} in the internal model.
    use, intrinsic :: ISO_C_Binding
    use FGSL
    use Numerical_Interpolation
    implicit none
    class           (intergalacticMediumStateInternal), intent(inout)   :: self
    double precision                                  , intent(in   )   :: time
    double precision                                  , dimension(0:1)  :: h
    integer         (c_size_t                        )                  :: i                              , j
    logical                                           , save            :: interpolationReset      =.true.
    type            (fgsl_interp_accel               ), save            :: interpolationAccelerator 
    !$omp threadprivate(interpolationReset,interpolationAccelerator)

    if (size(self%time) > 1) then
       i=Interpolate_Locate(self%time,interpolationAccelerator,time,reset=interpolationReset)
    else 
       internalFilteringMass=self%massFiltering(1)
       return
    end if
    if (self%time(i+1)-self%time(i) > 0.0d0) then
       h   =Interpolate_Linear_Generate_Factors(self%time,i,time)
    else
       h(0)=0.0d0 
       h(1)=1.0d0
    end if
    internalFilteringMass=0.0d0
    do j=0,1
       internalFilteringMass=internalFilteringMass+h(j)*self%massFiltering(i+j)
    end do
    return
  end function internalFilteringMass

  double precision function internalElectronFraction(self,time)
    !% Return the electron fraction of the \gls{igm} in the internal model.
    use, intrinsic :: ISO_C_Binding
    use FGSL
    use Numerical_Interpolation
    implicit none
    class           (intergalacticMediumStateInternal), intent(inout)   :: self
    double precision                                  , intent(in   )   :: time
    double precision                                  , dimension(0:1)  :: h
    integer         (c_size_t                        )                  :: i                              , j
    logical                                           , save            :: interpolationReset      =.true.
    type            (fgsl_interp_accel               ), save            :: interpolationAccelerator 
    !$omp threadprivate(interpolationReset,interpolationAccelerator)
    double precision                                                    :: densityHydrogen                , densityHelium

    if (size(self%time) > 1) then
       i=Interpolate_Locate(self%time,interpolationAccelerator,time,reset=interpolationReset)
    else 
       densityHydrogen         =self%densityHydrogen1(1)+self%densityHydrogen2(1)
       densityHelium           =self%densityHelium1  (1)+self%densityHelium2  (1)+self%densityHelium3(1)
       internalElectronFraction=+(self%densityHydrogen2(1)+self%densityHelium2(1)+2.0d0*self%densityHelium3(1)) &
            &                   /(densityHydrogen+2.0d0*densityHelium)
       return
    end if
    if (self%time(i+1)-self%time(i) > 0.0d0) then
       h   =Interpolate_Linear_Generate_Factors(self%time,i,time)
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
    !% Return the neutral hydrogen fraction of the \gls{igm} in the internal model.
    use, intrinsic :: ISO_C_Binding
    use FGSL
    use Numerical_Interpolation
    implicit none
    class           (intergalacticMediumStateInternal), intent(inout)   :: self
    double precision                                  , intent(in   )   :: time
    double precision                                  , dimension(0:1)  :: h
    logical                                           , save            :: interpolationReset      =.true.
    type            (fgsl_interp_accel               ), save            :: interpolationAccelerator 
    !$omp threadprivate(interpolationReset,interpolationAccelerator)
    integer         (c_size_t                        )                  :: i                              , j
    double precision                                                    :: densityHydrogen
    
    if (size(self%time) > 1) then
       i=Interpolate_Locate(self%time,interpolationAccelerator,time,reset=interpolationReset)
    else 
       densityHydrogen                =self%densityHydrogen1(1)+self%densityHydrogen2(1)
       internalNeutralHydrogenFraction=self%densityHydrogen1(1)/densityHydrogen
       return
    end if
    if (self%time(i+1)-self%time(i) > 0.0d0) then
       h   =Interpolate_Linear_Generate_Factors(self%time,i,time)
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
    !% Return the neutral helium fraction of the \gls{igm} in the internal model.
    use, intrinsic :: ISO_C_Binding
    use FGSL
    use Numerical_Interpolation
    implicit none
    class           (intergalacticMediumStateInternal), intent(inout)  :: self
    double precision                                  , intent(in   )  :: time
    double precision                                  , dimension(0:1) :: h
    logical                                           , save           :: interpolationReset      =.true.
    type            (fgsl_interp_accel               ), save           :: interpolationAccelerator 
    !$omp threadprivate(interpolationReset,interpolationAccelerator)
    integer         (c_size_t                        )                  :: i                              , j
    double precision                                                   :: densityHelium
    
    if (size(self%time) > 1) then
       i=Interpolate_Locate(self%time,interpolationAccelerator,time,reset=interpolationReset)
    else 
       densityHelium                =self%densityHelium1(1)+self%densityHelium2(1)+self%densityHelium3(1)
       internalNeutralHeliumFraction=self%densityHelium1(1)/densityHelium
       return
    end if
    if (self%time(i+1)-self%time(i) > 0.0d0) then
       h   =Interpolate_Linear_Generate_Factors(self%time,i,time)
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
    !% Return the singly ionized helium fraction of the \gls{igm} in the internal model.
    use, intrinsic :: ISO_C_Binding
    use FGSL
    use Numerical_Interpolation
    implicit none
    class           (intergalacticMediumStateInternal), intent(inout)   :: self
    double precision                                  , intent(in   )   :: time
    double precision                                  , dimension(0:1)  :: h
    logical                                           , save            :: interpolationReset      =.true.
    type            (fgsl_interp_accel)               , save            :: interpolationAccelerator 
    !$omp threadprivate(interpolationReset,interpolationAccelerator)
    integer         (c_size_t                        )                  :: i                              , j
    double precision                                                    :: densityHelium

    if (size(self%time) > 1) then
       i=Interpolate_Locate(self%time,interpolationAccelerator,time,reset=interpolationReset)
    else 
       densityHelium                      =self%densityHelium1(1)+self%densityHelium2(1)+self%densityHelium3(1)
       internalSinglyIonizedHeliumFraction=self%densityHelium2(1)/densityHelium
       return
    end if
    if (self%time(i+1)-self%time(i) > 0.0d0) then
       h   =Interpolate_Linear_Generate_Factors(self%time,i,time)
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
    !% Return the temperature of the \gls{igm} in the internal model.
    use, intrinsic :: ISO_C_Binding
    use FGSL
    use Numerical_Interpolation
    implicit none
    class           (intergalacticMediumStateInternal), intent(inout)   :: self
    double precision                                  , intent(in   )   :: time
    double precision                                  , dimension(0:1)  :: h
    logical                                           , save            :: interpolationReset      =.true.
    type            (fgsl_interp_accel)               , save            :: interpolationAccelerator 
    !$omp threadprivate(interpolationReset,interpolationAccelerator)
    integer         (c_size_t                        )                  :: i                              , j

    if (size(self%time) > 1) then
       i=Interpolate_Locate(self%time,interpolationAccelerator,time,reset=interpolationReset)
    else 
       internalTemperature=self%temperatureIGM(1)
       return
    end if
    if (self%time(i+1)-self%time(i) > 0.0d0) then
       h   =Interpolate_Linear_Generate_Factors(self%time,i,time)
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

  subroutine internalTimeSet(self,times)
    !% Set times in the internal intergalatic medium state class.
    use Memory_Management
    implicit none
    class           (intergalacticMediumStateInternal), intent(inout)               :: self
    double precision                                  , intent(in   ), dimension(:) :: times

    if (allocated(self%time)) call Dealloc_Array(self%time)
    call Alloc_Array(self%time,shape(times))
    self%time=times
    return
  end subroutine internalTimeSet

  subroutine internalTemperatureSet(self,temperatures)
    !% Set temperatures in the internal intergalatic medium state class.
    use Memory_Management
    implicit none
    class           (intergalacticMediumStateInternal), intent(inout)               :: self
    double precision                                  , intent(in   ), dimension(:) :: temperatures

    if (allocated(self%temperatureIGM)) call Dealloc_Array(self%temperatureIGM)
    call Alloc_Array(self%temperatureIGM,shape(temperatures))
    self%temperatureIGM=temperatures
    return
  end subroutine internalTemperatureSet

  subroutine internalDensityHydrogen1Set(self, densityHydrogen1)
    !% Set H1 densities in the internal intergalatic medium state class.
    use Memory_Management
    implicit none
    class           (intergalacticMediumStateInternal), intent(inout)               :: self
    double precision                                  , intent(in   ), dimension(:) :: densityHydrogen1

    if (allocated(self%densityHydrogen1)) call Dealloc_Array(self%densityHydrogen1)
    call Alloc_Array(self%densityHydrogen1,shape(densityHydrogen1))
    self%densityHydrogen1=densityHydrogen1
    return
  end subroutine internalDensityHydrogen1Set

  subroutine internalDensityHydrogen2Set(self, densityHydrogen2)
    !% Set H2 densities in the internal intergalatic medium state class.
    use Memory_Management
    implicit none
    class           (intergalacticMediumStateInternal), intent(inout)               :: self
    double precision                                  , intent(in   ), dimension(:) :: densityHydrogen2

    if (allocated(self%densityHydrogen2)) call Dealloc_Array(self%densityHydrogen2)
    call Alloc_Array(self%densityHydrogen2,shape(densityHydrogen2))
    self%densityHydrogen2=densityHydrogen2
    return
  end subroutine internalDensityHydrogen2Set

  subroutine internalDensityHelium1Set(self,densityHelium1)
    !% Set He1 densities in the internal intergalatic medium state class.
    use Memory_Management
    implicit none
    class           (intergalacticMediumStateInternal), intent(inout)               :: self
    double precision                                  , intent(in   ), dimension(:) :: densityHelium1

    if (allocated(self%densityHelium1)) call Dealloc_Array(self%densityHelium1)
    call Alloc_Array(self%densityHelium1,shape(densityHelium1))
    self%densityHelium1=densityHelium1
    return
  end subroutine internalDensityHelium1Set

  subroutine internalDensityHelium2Set(self,densityHelium2)
    !% Set He2 densities in the internal intergalatic medium state class.
    use Memory_Management
    implicit none
    class           (intergalacticMediumStateInternal), intent(inout)               :: self
    double precision                                  , intent(in   ), dimension(:) :: densityHelium2

    if (allocated(self%densityHelium2)) call Dealloc_Array(self%densityHelium2)
    call Alloc_Array(self%densityHelium2,shape(densityHelium2))
    self%densityHelium2=densityHelium2
    return
  end subroutine internalDensityHelium2Set

  subroutine internalDensityHelium3Set(self, densityHelium3)
    !% Set He3 densities in the internal intergalatic medium state class.
    use Memory_Management
    implicit none
    class           (intergalacticMediumStateInternal), intent(inout)               :: self
    double precision                                  , intent(in   ), dimension(:) :: densityHelium3

    if (allocated(self%densityHelium3)) call Dealloc_Array(self%densityHelium3)
    call Alloc_Array(self%densityHelium3,shape(densityHelium3))
    self%densityHelium3=densityHelium3
    return
  end subroutine internalDensityHelium3Set

  subroutine internalMassFilteringSet(self,massFiltering)
    !% Set filtering masses in the internal intergalatic medium state class.
    use Memory_Management
    implicit none
    class           (intergalacticMediumStateInternal), intent(inout)               :: self
    double precision                                  , intent(in   ), dimension(:) :: massFiltering

    if (allocated(self%massFiltering)) call Dealloc_Array(self%massFiltering)
    call Alloc_Array(self%massFiltering,shape(massFiltering))
    self%massFiltering=massFiltering
    return
  end subroutine internalMassFilteringSet

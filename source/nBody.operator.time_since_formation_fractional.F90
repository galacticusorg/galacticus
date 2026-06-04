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
Implements an N-body data operator which computes the time since formation as a fraction of crossing time for halos.
!!}
  
  use :: Cosmology_Functions, only : cosmologyFunctionsClass

  !![
  <nbodyOperator name="nbodyOperatorTimeSinceFormationFractional">
   <description>An N-body data operator which computes and stores the time since formation as a fraction of crossing time for halos.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorTimeSinceFormationFractional
     !!{
     An N-body data operator which computes the time since formation as a fraction of crossing time for halos.
     !!}
     private
     class(cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null()
   contains
     final     ::            timeSinceFormationDestructor
     procedure :: operate => timeSinceFormationFractionalOperate
  end type nbodyOperatorTimeSinceFormationFractional

  interface nbodyOperatorTimeSinceFormationFractional
     !!{
     Constructors for the \refClass{nbodyOperatorTimeSinceFormationFractional} N-body operator class.
     !!}
     module procedure timeSinceFormationFractionalConstructorParameters
     module procedure timeSinceFormationConstructorInternal
  end interface nbodyOperatorTimeSinceFormationFractional

contains

  function timeSinceFormationFractionalConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{nbodyOperatorTimeSinceFormationFractional} N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nbodyOperatorTimeSinceFormationFractional)                :: self
    type (inputParameters                          ), intent(inout) :: parameters
    class(cosmologyFunctionsClass                  ), pointer       :: cosmologyFunctions_

    !![
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    self=nbodyOperatorTimeSinceFormationFractional(cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function timeSinceFormationFractionalConstructorParameters

  function timeSinceFormationConstructorInternal(cosmologyFunctions_) result (self)
    !!{
    Internal constructor for the ``timeSinceFormation'' N-body operator class.
    !!}
    implicit none
    type (nbodyOperatorTimeSinceFormationFractional)                        :: self
    class(cosmologyFunctionsClass                  ), intent(in   ), target :: cosmologyFunctions_
    !![
    <constructorAssign variables="*cosmologyFunctions_"/>
    !!]

    return
  end function timeSinceFormationConstructorInternal
  
  subroutine timeSinceFormationDestructor(self)
    !!{
    Destructor for the ``timeSinceFormation'' N-body operator class.
    !!}
    implicit none
    type(nbodyOperatorTimeSinceFormationFractional), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine timeSinceFormationDestructor

  subroutine timeSinceFormationFractionalOperate(self,simulations)
    !!{
    Compute the distance of each particle from a point.
    !!}
    use    :: Display      , only : displayIndent      , displayUnindent, verbosityLevelStandard, displayCounter, &
         &                          displayCounterClear
    !$ use :: OMP_Lib      , only : OMP_Get_Thread_Num
#ifdef USEMPI
    use    :: MPI_Utilities, only : mpiSelf
#endif
    implicit none
    class           (nbodyOperatorTimeSinceFormationFractional), intent(inout)                 :: self
    type            (nBodyData                                ), intent(inout), dimension(:  ) :: simulations
    class           (cosmologyFunctionsClass                  ), pointer                       :: cosmologyFunctions_
    double precision                                           , pointer      , dimension(:  ) :: expansionFactor    , timeFormation               , &
         &                                                                                        timeCrossing       , timeSinceFormationFractional
    integer         (c_size_t                                 )                                :: iSimulation        , i

    call displayIndent('compute times since halo formation',verbosityLevelStandard)
    do iSimulation=1,size(simulations)
       ! Retrieve required properties.
       expansionFactor => simulations(iSimulation)%propertiesReal%value('expansionFactor')
       timeFormation   => simulations(iSimulation)%propertiesReal%value('timeFormation'  )
       timeCrossing    => simulations(iSimulation)%propertiesReal%value('timeCrossing'   )
       ! Allocate workspace.
       allocate(timeSinceFormationFractional(size(timeFormation)))
       ! Compute times since formation in units of the crossing time.
       !$omp parallel private(i,cosmologyFunctions_)
       allocate(cosmologyFunctions_,mold=self%cosmologyFunctions_)
       !$omp critical(nBodyOperatorTimeFormationDeepCopy)
       !![
       <deepCopyReset variables="self%cosmologyFunctions_"/>
       <deepCopy source="self%cosmologyFunctions_" destination="cosmologyFunctions_"/>
       <deepCopyFinalize variables="cosmologyFunctions_"/>
       !!]
       !$omp end critical(nBodyOperatorTimeFormationDeepCopy)
       !$omp barrier
       ! Visit all particles.
       !$omp do schedule(dynamic)
       do i=1_c_size_t,size(expansionFactor)
          timeSinceFormationFractional(i)=+(                                                       &
               &                            +cosmologyFunctions_%cosmicTime   (expansionFactor(i)) &
               &                            -                    timeFormation                (i)  &
               &                           )                                                       &
               &                          /                      timeCrossing                 (i)
          ! Update progress.
          !$ if (OMP_Get_Thread_Num() == 0) then
#ifdef USEMPI
          if (mpiSelf%isMaster()) then
#endif
             call displayCounter(                                 &
                  &              int(                             &
                  &                  +100.0d0                     &
                  &                  *float(i                  )  &
                  &                  /float(size(timeFormation))  &
                  &                 )                           , &
                  &              .false.                          &
                  &             )
#ifdef USEMPI
          end if
#endif
          !$ end if
       end do
       !$omp end do
       deallocate(cosmologyFunctions_)
       !$omp end parallel
#ifdef USEMPI
       if (mpiSelf%isMaster()) then
#endif
          call displayCounterClear()
#ifdef USEMPI
       end if
#endif
       ! Store results.
       call simulations(iSimulation)%propertiesReal%set('timeSinceFormationFractional',timeSinceFormationFractional)
       nullify(timeSinceFormationFractional)
    end do
    call displayUnindent('done',verbosityLevelStandard)
    return
  end subroutine timeSinceFormationFractionalOperate

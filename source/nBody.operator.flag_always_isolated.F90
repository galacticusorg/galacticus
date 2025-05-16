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

!!{
Implements an N-body data operator which flags particles that have been always isolated.
!!}

  !![
  <nbodyOperator name="nbodyOperatorFlagAlwaysIsolated">
   <description>An N-body data operator which flags particles that have been always isolated.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorFlagAlwaysIsolated
     !!{
     An N-body data operator which flags particles that have been always isolated.
     !!}
     private
     double precision :: massFactor
   contains
     procedure :: operate => flagAlwaysIsolatedOperate
  end type nbodyOperatorFlagAlwaysIsolated

  interface nbodyOperatorFlagAlwaysIsolated
     !!{
     Constructors for the \refClass{nbodyOperatorFlagAlwaysIsolated} N-body operator class.
     !!}
     module procedure flagAlwaysIsolatedConstructorParameters
     module procedure flagAlwaysIsolatedConstructorInternal
  end interface nbodyOperatorFlagAlwaysIsolated

contains

  function flagAlwaysIsolatedConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{nbodyOperatorFlagAlwaysIsolated} N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nbodyOperatorFlagAlwaysIsolated)                :: self
    type            (inputParameters                ), intent(inout) :: parameters
    double precision                                                 :: massFactor

    !![
    <inputParameter>
      <name>massFactor</name>
      <source>parameters</source>
      <description>The factor by which virial mass must increase for previous non-isolation to be ignored.</description>
      <defaultValue>2.0d0</defaultValue>
    </inputParameter>
    !!]
    self=nbodyOperatorFlagAlwaysIsolated(massFactor)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function flagAlwaysIsolatedConstructorParameters

  function flagAlwaysIsolatedConstructorInternal(massFactor) result (self)
    !!{
    Internal constructor for the \refClass{nbodyOperatorFlagAlwaysIsolated} N-body operator class.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nbodyOperatorFlagAlwaysIsolated)                :: self
    double precision                                 , intent(in   ) :: massFactor
    !![
    <constructorAssign variables="massFactor"/>
    !!]

    return
  end function flagAlwaysIsolatedConstructorInternal

  subroutine flagAlwaysIsolatedOperate(self,simulations)
    !!{
    Identify and flag particles which have been always isolated.
    !!}
    use    :: Arrays_Search     , only : searchIndexed
    use    :: Display           , only : displayCounter         , displayCounterClear, displayIndent, displayUnindent, &
          &                              verbosityLevelStandard
    use    :: Error             , only : Error_Report
    use    :: ISO_Varying_String, only : var_str
    !$ use :: OMP_Lib           , only : OMP_Get_Thread_Num
    use    :: Sorting           , only : sortIndex
    use    :: String_Handling   , only : operator(//)
    implicit none
    class           (nbodyOperatorFlagAlwaysIsolated), intent(inout)               :: self
    type            (nBodyData                      ), intent(inout), dimension(:) :: simulations
    integer         (c_size_t                       ), pointer      , dimension(:) :: alwaysIsolated         , hostID     , &
         &                                                                            descendantID           , particleIDs, &
         &                                                                            isMostMassiveProgenitor
    integer         (c_size_t                       ), allocatable  , dimension(:) :: indexID
    double precision                                 , pointer      , dimension(:) :: massVirial
    double precision                                                               :: massReference
    integer         (c_size_t                       )                              :: i                      , j          , &
         &                                                                            k                      , iSimulation

    call displayIndent('flag always isolated objects',verbosityLevelStandard)
    do iSimulation=1,size(simulations)
       ! Retrieve required properties.
       particleIDs             => simulations(iSimulation)%propertiesInteger%value('particleID'             )
       hostID                  => simulations(iSimulation)%propertiesInteger%value('hostID'                 )
       descendantID            => simulations(iSimulation)%propertiesInteger%value('descendantID'           )
       isMostMassiveProgenitor => simulations(iSimulation)%propertiesInteger%value('isMostMassiveProgenitor')
       massVirial              => simulations(iSimulation)%propertiesReal   %value('massVirial'             )
       ! Allocate workspace.
       allocate(indexID       (size(particleIDs)))
       allocate(alwaysIsolated(size(particleIDs)))
       ! Build a sort index.
       indexID=sortIndex(particleIDs)
       ! Initialize status - assuming all particles are always-isolated initially.
       alwaysIsolated=1_c_size_t
       ! Visit each particle.
       !$omp parallel do private(j,k,massReference) schedule(dynamic)
       do i=1_c_size_t,size(alwaysIsolated)
          ! Skip isolated halos.
          if (hostID(i) < 0_c_size_t) cycle
          ! Trace descendants, marking as not-always-isolated until mass increases sufficiently.
          j            =           i
          massReference=massVirial(i)
          do while (massVirial(j) < self%massFactor*massReference)             
             alwaysIsolated(j)=0_c_size_t
             if (descendantID(j) < 0_c_size_t .or. isMostMassiveProgenitor(j) == 0) exit
             k=searchIndexed(particleIDs,indexID,descendantID(j))
             if     (                                           &
                  &   k                         < 1_c_size_t    &
                  &  .or.                                       &
                  &   k                         > size(indexID) &
                  & )                                           &
                  & call Error_Report('failed to find descendant'//{introspection:location})
             k=indexID(k)
             if     (                                   &
                  &   particleIDs(k) /= descendantID(j) &
                  & )                                   &
                  & call Error_Report(var_str('failed to find descendant [')//descendantID(j)//'] of ['//particleIDs(j)//']'//{introspection:location})
             j=k
          end do
          !$ if (OMP_Get_Thread_Num() == 0) then
          call displayCounter(int(100.0d0*dble(i)/dble(size(alwaysIsolated))),verbosity=verbosityLevelStandard,isNew=i == 1_c_size_t)
          !$ end if
       end do
       !$omp end parallel do
       call displayCounterClear(verbosityLevelStandard)
       ! Store results.
       call simulations(iSimulation)%propertiesInteger%set('alwaysIsolated',alwaysIsolated)
       deallocate(indexID       )
    end do
    call displayUnindent('done',verbosityLevelStandard)
    return
  end subroutine flagAlwaysIsolatedOperate

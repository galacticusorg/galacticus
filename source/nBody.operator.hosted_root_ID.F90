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
  Implements an N-body data operator which determines an ID of the root halo found by following hosts.
  !!}

  !![
  <nbodyOperator name="nbodyOperatorHostedRootID">
   <description>An N-body data operator which determines an ID of the root halo found by following hosts.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorHostedRootID
     !!{
     An N-body data operator which determines an ID of the root halo found by following hosts.
     !!}
     private
     logical :: missingHalosAreFatal
   contains
     procedure :: operate => hostedRootIDOperate
  end type nbodyOperatorHostedRootID

  interface nbodyOperatorHostedRootID
     !!{
     Constructors for the \refClass{nbodyOperatorHostedRootID} N-body operator class.
     !!}
     module procedure hostedRootIDConstructorParameters
     module procedure hostedRootIDConstructorInternal
  end interface nbodyOperatorHostedRootID

contains

  function hostedRootIDConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{nbodyOperatorHostedRootID} N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (nbodyOperatorHostedRootID)                :: self
    type   (inputParameters          ), intent(inout) :: parameters
    logical                                           :: missingHalosAreFatal

    !![
    <inputParameter>
      <name>missingHalosAreFatal</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>If true, if a halo is not found during the search through hosts and descendants then a fatal error occurs. Otherwise, such missing halos are ignored, and a {\normalfont \ttfamily hostedRootID} value of $-1$ is assigned to the particle.</description>
    </inputParameter>
    !!]
    self=nbodyOperatorHostedRootID(missingHalosAreFatal)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function hostedRootIDConstructorParameters

  function hostedRootIDConstructorInternal(missingHalosAreFatal) result (self)
    !!{
    Internal constructor for the \refClass{nbodyOperatorHostedRootID} N-body operator class.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (nbodyOperatorHostedRootID)                :: self
    logical                           , intent(in   ) :: missingHalosAreFatal
    !![
    <constructorAssign variables="missingHalosAreFatal"/>
    !!]
    
    return
  end function hostedRootIDConstructorInternal

  subroutine hostedRootIDOperate(self,simulations)
    !!{
    Determine an ID of the root halo found by following hosts.
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
    class  (nbodyOperatorHostedRootID), intent(inout)               :: self
    type   (nBodyData                ), intent(inout), dimension(:) :: simulations
    integer(c_size_t                 ), pointer      , dimension(:) :: hostedRootID, isolatedHostID, &
         &                                                             descendantID, particleIDs
    integer(c_size_t                 ), allocatable  , dimension(:) :: indexID
    integer(c_size_t                 )                              :: i           , j             , &
         &                                                             k           , l             , &
         &                                                             iSimulation

    call displayIndent('determine hosted root IDs',verbosityLevelStandard)
    do iSimulation=1,size(simulations)
       ! Retrieve required properties.
       particleIDs    => simulations(iSimulation)%propertiesInteger%value('particleID'    )
       isolatedHostID => simulations(iSimulation)%propertiesInteger%value('isolatedHostID')
       descendantID   => simulations(iSimulation)%propertiesInteger%value('descendantID'  )
       ! Allocate workspace.
       allocate(indexID     (size(particleIDs)))
       allocate(hostedRootID(size(particleIDs)))
       ! Build a sort index.
       indexID=sortIndex(particleIDs)
       ! Visit each particle.
       !$omp parallel do private(j,k,l) schedule(dynamic)
       do i=1_c_size_t,size(hostedRootID)
          ! Trace descendants through hosts.
          !! Begin at the current particle.
          j=i
          !! Step up through hosts and down through descendants until none remain.
          do while (descendantID(j) >= 0_c_size_t .or. isolatedHostID(j) >= 0_c_size_t)
             !! Move up through hosts.
             if      (isolatedHostID(j) >= 0_c_size_t) then
                l=isolatedHostID(j)
             else if (descendantID  (j) > 0_c_size_t) then
                l=descendantID  (j)
             else
                call Error_Report('no host or descendant - this should not happen'//{introspection:location})
             end if
             k=searchIndexed(particleIDs,indexID,l)
             if     (                                 &
                  &   k               < 1_c_size_t    &
                  &  .or.                             &
                  &   k               > size(indexID) &
                  & )                                 &
                  & then
                if (self%missingHalosAreFatal) then
                   call Error_Report('failed to find next halo'//{introspection:location})
                else
                   j=-1_c_size_t
                   exit
                end if
             end if
             k=indexID(k)
             if     (                                 &
                  &   particleIDs(k) /= l             &
                  & )                                 &
                  & then
                if (self%missingHalosAreFatal) then
                   call Error_Report(var_str('failed to find next halos [')//l//'] of ['//particleIDs(j)//']'//{introspection:location})
                else
                   j=-1_c_size_t
                   exit
                end if
             end if
             j=k
          end do
          ! Assign the hosted root ID.
          if (j < 0_c_size_t) then
             ! No hosted root was found (due to a missing halo). Assign a negative value to indicate this.
             hostedrootid(i)=-1
          else
             ! A hosted root was found, record its ID.
             hostedRootID(i)=particleIDs(j)
          end if
          !$ if (OMP_Get_Thread_Num() == 0) then
          call displayCounter(int(100.0d0*dble(i)/dble(size(hostedRootID))),verbosity=verbosityLevelStandard,isNew=i == 1_c_size_t)
          !$ end if
       end do
       !$omp end parallel do
       call displayCounterClear(verbosityLevelStandard)
       ! Store results.
       call simulations(iSimulation)%propertiesInteger%set('hostedRootID',hostedRootID)
       deallocate(indexID     )
       nullify   (hostedRootID)
    end do
    call displayUnindent('done',verbosityLevelStandard)
    return
  end subroutine hostedRootIDOperate

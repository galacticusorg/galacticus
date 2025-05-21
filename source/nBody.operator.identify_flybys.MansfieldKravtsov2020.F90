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
Implements an N-body data operator which identifies flyby halos following the algorithm of \cite{mansfield_three_2020}.
!!}

  !![
  <nbodyOperator name="nbodyOperatorIdentifyFlybysMansfieldKravtsov2020">
   <description>An N-body data operator which identifies flyby halos following the algorithm of \cite{mansfield_three_2020}.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorIdentifyFlybysMansfieldKravtsov2020
     !!{
     An N-body data operator which identifies flyby halos following the algorithm of \cite{mansfield_three_2020}.
     !!}
     private
     logical :: missingHostsAreFatal
   contains
     procedure :: operate => identifyFlybysMansfieldKravtsov2020Operate
  end type nbodyOperatorIdentifyFlybysMansfieldKravtsov2020

  interface nbodyOperatorIdentifyFlybysMansfieldKravtsov2020
     !!{
     Constructors for the \refClass{nbodyOperatorIdentifyFlybysMansfieldKravtsov2020} N-body operator class.
     !!}
     module procedure identifyFlybysMansfieldKravtsov2020ConstructorParameters
     module procedure identifyFlybysMansfieldKravtsov2020ConstructorInternal
  end interface nbodyOperatorIdentifyFlybysMansfieldKravtsov2020

contains

  function identifyFlybysMansfieldKravtsov2020ConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{nbodyOperatorIdentifyFlybysMansfieldKravtsov2020} N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (nbodyOperatorIdentifyFlybysMansfieldKravtsov2020)                :: self
    type   (inputParameters                                 ), intent(inout) :: parameters
    logical                                                                  :: missingHostsAreFatal
    
    !![
    <inputParameter>
      <name>missingHostsAreFatal</name>
      <source>parameters</source>
      <description>If true, missing hosts cause a fatal error. Otherwise, a missing host causes the halo to be ignored in flyby analysis of its descendants.</description>
    </inputParameter>
    !!]
    self=nbodyOperatorIdentifyFlybysMansfieldKravtsov2020(missingHostsAreFatal)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function identifyFlybysMansfieldKravtsov2020ConstructorParameters

  function identifyFlybysMansfieldKravtsov2020ConstructorInternal(missingHostsAreFatal) result (self)
    !!{
    Internal constructor for the \refClass{nbodyOperatorIdentifyFlybysMansfieldKravtsov2020} N-body operator class.
    !!}
    implicit none
    type   (nbodyOperatorIdentifyFlybysMansfieldKravtsov2020)                :: self
    logical                                                  , intent(in   ) :: missingHostsAreFatal
    !![
    <constructorAssign variables="missingHostsAreFatal"/>
    !!]
    
    return
  end function identifyFlybysMansfieldKravtsov2020ConstructorInternal

  subroutine identifyFlybysMansfieldKravtsov2020Operate(self,simulations)
    !!{
    Identify flyby halos.
    !!}
    use    :: Arrays_Search     , only : searchIndexed
    use    :: Display           , only : displayCounter         , displayCounterClear, displayIndent, displayUnindent, &
          &                              verbosityLevelStandard
    use    :: Error             , only : Error_Report
    use    :: ISO_Varying_String, only : var_str
    !$ use :: OMP_Lib           , only : OMP_Get_Thread_Num
    use    :: Sorting           , only : sortIndex
    use    :: String_Handling   , only : operator(//)
    use    :: Vectors           , only : Vector_Magnitude
    implicit none
    class           (nbodyOperatorIdentifyFlybysMansfieldKravtsov2020), intent(inout)                 :: self
    type            (nBodyData                                       ), intent(inout), dimension(:  ) :: simulations
    integer         (c_size_t                                        ), pointer      , dimension(:  ) :: isFlyby                , hostID      , &
         &                                                                                               descendantID           , particleIDs , &
         &                                                                                               isMostMassiveProgenitor, snapshotID
    integer         (c_size_t                                        ), allocatable  , dimension(:  ) :: indexID
    double precision                                                  , pointer      , dimension(:  ) :: massVirial             , radiusVirial, &
         &                                                                                               expansionFactor
    double precision                                                  , pointer      , dimension(:,:) :: position
    double precision                                                                 , dimension(3  ) :: separation
    integer         (c_size_t                                        )                                :: i                      , iSimulation , &
         &                                                                                               jHalo                  , jHost       , &
         &                                                                                               kHalo                  , kHost       , &
         &                                                                                               iSimulation            , m
    double precision                                                                                  :: boxSize
    logical                                                                                           :: hostDescendantExists
    
    call displayIndent('identify flyby halos',verbosityLevelStandard)
    do iSimulation=1,size(simulations)
       ! Get box size.
       if (simulations(iSimulation)%attributesReal%exists('boxSize')) then
          boxSize=simulations(iSimulation)%attributesReal%value('boxSize')
       else
          boxSize=-1.0d0
       end if
       ! Retrieve required properties.
       particleIDs             => simulations(iSimulation)%propertiesInteger  %value('particleID'             )
       hostID                  => simulations(iSimulation)%propertiesInteger  %value('hostID'                 )
       descendantID            => simulations(iSimulation)%propertiesInteger  %value('descendantID'           )
       isMostMassiveProgenitor => simulations(iSimulation)%propertiesInteger  %value('isMostMassiveProgenitor')
       snapshotID              => simulations(iSimulation)%propertiesInteger  %value('snapshotID'             )  
       radiusVirial            => simulations(iSimulation)%propertiesReal     %value('radiusVirial'           )
       massVirial              => simulations(iSimulation)%propertiesReal     %value('massVirial'             )
       expansionFactor         => simulations(iSimulation)%propertiesReal     %value('expansionFactor'        )
       position                => simulations(iSimulation)%propertiesRealRank1%value('position'               )
       ! Allocate workspace.
       allocate(indexID(size(particleIDs)))
       allocate(isFlyby(size(particleIDs)))
       ! Build a sort index.
       indexID=sortIndex(particleIDs)
       ! Initialize status - assuming all particles are always-isolated initially.
       isFlyby=0_c_size_t
       ! Visit each particle.
       !$omp parallel do private(jHost,jHalo,kHost,kHalo,hostDescendantExists,m,separation) schedule(dynamic)
       do i=1_c_size_t,size(isFlyby)
          ! Skip isolated halos.
          if (hostID(i) < 0_c_size_t) cycle
          ! Find the host halo.
          kHost=searchIndexed(particleIDs,indexID,hostID(i))
          if     (                       &
               &   kHost < 1_c_size_t    &
               &  .or.                   &
               &   kHost > size(indexID) &
               & )                       &
               & then
             if (self%missingHostsAreFatal) then
                call Error_Report('failed to find host'//{introspection:location})
             else
                cycle
             end if
          end if
          kHost=indexID(kHost)
          if     (                                 &
               &   particleIDs(kHost) /= hostID(i) &
               & )                                 &
               & then
             if (self%missingHostsAreFatal) then
                call Error_Report(var_str('failed to find host [')//hostID(i)//'] of ['//particleIDs(i)//']'//{introspection:location})
             else
                cycle
             end if
          end if
          jHost=kHost
          ! Trace descendants, marking as flybys.
          jHalo               =i
          hostDescendantExists=.true.
          do while (hostDescendantExists)
             isFlyby(jHalo)=1_c_size_t
             ! If there is no descendant, or the current halo is not the most massive progenitor, we are done marking flybys, so exit.
             if (descendantID(jHalo) < 0_c_size_t .or. isMostMassiveProgenitor(jHalo) == 0) exit
             ! Find the descendant halo.
             kHalo=searchIndexed(particleIDs,indexID,descendantID(jHalo))
             if     (                       &
                  &   kHalo < 1_c_size_t    &
                  &  .or.                   &
                  &   kHalo > size(indexID) &
                  & )                       &
                  & call Error_Report('failed to find descendant'//{introspection:location})
             kHalo=indexID(kHalo)
             if     (                                           &
                  &   particleIDs(kHalo) /= descendantID(jHalo) &
                  & )                                           &
                  & call Error_Report(var_str('failed to find descendant [')//descendantID(jHalo)//'] of ['//particleIDs(jHalo)//']'//{introspection:location})
             jHalo=kHalo
             ! We must now find the descendant of the host halo at this same snapshot.
             do while (snapshotID(jHost) < snapshotID(jHalo))
                if (descendantID(jHost) < 0_c_size_t) then
                   ! No descendant. Record this and exit.
                   hostDescendantExists=.false.
                   exit
                end if
                kHost=searchIndexed(particleIDs,indexID,descendantID(jHost))
                if     (                       &
                     &   kHost < 1_c_size_t    &
                     &  .or.                   &
                     &   kHost > size(indexID) &
                     & )                       &
                     & call Error_Report('failed to find host'//{introspection:location})
                kHost=indexID(kHost)
                if     (                                           &
                     &   particleIDs(kHost) /= descendantID(jHost) &
                     & )                                           &
                     & call Error_Report(var_str('failed to find descendant [')//descendantID(jHost)//'] of ['//particleIDs(jHost)//']'//{introspection:location})
                jHost=kHost
                if (snapshotID(jHost) > snapshotID(jHalo)) then
                   ! We have jumped past the snapshot of the halo, so the host did not have a descendant at the relevant snapshot
                   ! - most likely because that descendant was sub-resolution. Record this and exit.
                   hostDescendantExists=.false.
                   exit
                end if
             end do
             ! Check conditions from Mansfield & Kravtsov (2020; §2.3.1).
             !! (  i) Host must have a descendant at our snapshot.
             if (.not.hostDescendantExists) exit
             !! ( ii) Descendant of the host must not be within the virial radius of descendant of our halo.
             separation=position(:,jHalo)-position(:,jHost)
             if (boxSize > 0.0d0) then
                do m=1,3
                   if (separation(m) > +boxSize*expansionFactor(jHalo)/2.0d0) separation(m)=separation(m)-boxSize*expansionFactor(jHalo)
                   if (separation(m) < -boxSize*expansionFactor(jHalo)/2.0d0) separation(m)=separation(m)+boxSize*expansionFactor(jHalo)
                end do                
             end if
             if (Vector_Magnitude(separation) < radiusVirial(jHalo)) exit
             !! (iii) Descendant of host must have higher mass than descendant of our halo.
             if (massVirial(jHost) <= massVirial(jHalo)) exit
          end do
          !$ if (OMP_Get_Thread_Num() == 0) then
          call displayCounter(int(100.0d0*dble(i)/dble(size(isFlyby))),verbosity=verbosityLevelStandard,isNew=i == 1_c_size_t)
          !$ end if
       end do
       !$omp end parallel do
       call displayCounterClear(verbosityLevelStandard)
       ! Store results.
       call simulations(iSimulation)%propertiesInteger%set('isFlyby',isFlyby)
       deallocate(indexID)
    end do
    call displayUnindent('done',verbosityLevelStandard)
    return
  end subroutine identifyFlybysMansfieldKravtsov2020Operate

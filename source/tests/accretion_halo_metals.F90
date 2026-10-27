!! Copyright 2026
!!    Andrew Robertson
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

!!{RST
Tests that initial halo metal masses scale with gas masses, independently of the instantaneous accretion rate.
!!}

!+ Contributions to this file made by: Andrew Robertson, ChatGPT.

program Test_Accretion_Halo_Metals
  !!{RST
  Check accreted and failed-to-accrete metal masses against the gas reservoir masses, and metal rates against the gas rates.
  !!}
  use :: Abundances_Structure            , only : abundances
  use :: Accretion_Halos                 , only : accretionHaloClass                 , accretionModeCold                , &
  &                                               accretionModeHot                   , accretionModeTotal               , &
  &                                               enumerationAccretionModeType
  use :: Display                         , only : displayVerbositySet                , verbosityLevelStandard
  use :: Error                           , only : Error_Handler_Register
  use :: Events_Hooks                    , only : eventsHooksInitialize
  use :: Functions_Global_Utilities      , only : Functions_Global_Set
  use :: Galacticus_Nodes                , only : nodeClassHierarchyFinalize         , nodeClassHierarchyInitialize     , &
  &                                               nodeComponentBasic                 , nodeComponentHotHalo             , &
  &                                               treeNode
  use :: Input_Parameters                , only : inputParameters
  use :: Node_Components                 , only : Node_Components_Initialize         , Node_Components_Thread_Initialize, &
  &                                               Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Numerical_Constants_Astronomical, only : metallicitySolar
  use :: Unit_Tests                      , only : Assert                             , Unit_Tests_Begin_Group           , &
  &                                               Unit_Tests_End_Group               , Unit_Tests_Finish
  implicit none
  type            (inputParameters             )               :: parameters
  class           (accretionHaloClass          ), pointer      :: accretionHalo_
  type            (treeNode                    ), pointer      :: node
  class           (nodeComponentBasic          ), pointer      :: basic
  class           (nodeComponentHotHalo        ), pointer      :: hotHalo
  type            (abundances                  )               :: metals
  type            (enumerationAccretionModeType), dimension(3) :: modes
  ! The fixed metallicity is 0.1 Solar and the baryon fraction is 0.045/0.3, as specified in the parameter file.
  double precision                              , parameter    :: metallicityIGM    =0.1d0*metallicitySolar, &
       &                                                          baryonFraction    =0.15d0                , &
       &                                                          massHalo          =1.0d10
  ! Before and after reionization. The deliberately large suppression velocity ensures that the latter case puts all gas
  ! into the failed reservoir. Both mass methods therefore have a non-zero expected result in at least one configuration.
  double precision                              , parameter    :: times(2)          =[0.1d0,10.0d0]        , &
       &                                                          rates(3)          =[3.0d9,0.0d0,-7.0d9]
  double precision                                             :: fractionFailed                           , &
       &                                                          massExpected                             , &
       &                                                          massFailedExpected                       , &
       &                                                          rateExpected                             , &
       &                                                          rateFailedExpected
  integer                                                      :: iTime                                    , &
       &                                                          iRate                                    , &
       &                                                          iMode
  character       (len=80                      )               :: label

  call displayVerbositySet  (verbosityLevelStandard)
  call Error_Handler_Register()
  call eventsHooksInitialize()
  call Functions_Global_Set ()
  parameters=inputParameters('testSuite/parameters/accretionHaloMetals.xml')
  call nodeClassHierarchyInitialize     (parameters)
  call Node_Components_Initialize       (parameters)
  call Node_Components_Thread_Initialize(parameters)
  !![
  <objectBuilder class="accretionHalo" name="accretionHalo_" source="parameters"/>
  !!]
  modes=[accretionModeHot,accretionModeTotal,accretionModeCold]
  call Unit_Tests_Begin_Group("Halo accretion: metal masses and rates")
  do iTime=1,size(times)
     fractionFailed=real(iTime-1,kind=kind(0.0d0))
     do iRate=1,size(rates)
        ! A fresh node avoids carrying cached halo properties from one time/rate configuration to another.
        node    => treeNode     (                 )
        basic   => node%basic   (autoCreate=.true.)
        hotHalo => node%hotHalo(autoCreate=.true.)
        call basic  %massSet            (massHalo    )
        call basic  %timeSet            (times(iTime))
        call basic  %timeLastIsolatedSet(times(iTime))
        call basic  %accretionRateSet   (rates(iRate))
        call hotHalo%unaccretedMassSet  (0.0d0       )
        do iMode=1,size(modes)
           write(label,'(a,i0,a,i0,a,i0)') 'time ',iTime,', rate ',iRate,', mode ',iMode
           massExpected      =baryonFraction*massHalo    *(1.0d0-fractionFailed)
           massFailedExpected=baryonFraction*massHalo    *       fractionFailed
           rateExpected      =baryonFraction*rates(iRate)*(1.0d0-fractionFailed)
           rateFailedExpected=baryonFraction*rates(iRate)*       fractionFailed
           if (iMode == 3) then
              massExpected      =0.0d0
              massFailedExpected=0.0d0
              rateExpected      =0.0d0
              rateFailedExpected=0.0d0
           end if
           ! Check the reservoir partition independently of the metal methods, including the empty cold mode.
           call Assert(trim(label)//': accreted gas mass',accretionHalo_%accretedMass      (node,modes(iMode)),massExpected      ,relTol=1.0d-12,absTol=1.0d-12)
           call Assert(trim(label)//': failed gas mass'  ,accretionHalo_%failedAccretedMass(node,modes(iMode)),massFailedExpected,relTol=1.0d-12,absTol=1.0d-12)
           ! These are masses, even when the halo is losing mass or its instantaneous accretion rate is zero.
           metals=accretionHalo_%accretedMassMetals(node,modes(iMode))
           call Assert(trim(label)//': accreted metal mass',metals%metallicity(),metallicityIGM*massExpected      ,relTol=1.0d-12,absTol=1.0d-12)
           metals=accretionHalo_%failedAccretedMassMetals(node,modes(iMode))
           call Assert(trim(label)//': failed metal mass'  ,metals%metallicity(),metallicityIGM*massFailedExpected,relTol=1.0d-12,absTol=1.0d-12)
           ! The rate methods must retain their rate dependence, including negative rates; the fix must not zero them.
           metals=accretionHalo_%accretionRateMetals(node,modes(iMode))
           call Assert(trim(label)//': accreted metal rate',metals%metallicity(),metallicityIGM*rateExpected      ,relTol=1.0d-12,absTol=1.0d-12)
           metals=accretionHalo_%failedAccretionRateMetals(node,modes(iMode))
           call Assert(trim(label)//': failed metal rate'  ,metals%metallicity(),metallicityIGM*rateFailedExpected,relTol=1.0d-12,absTol=1.0d-12)
        end do
        call node%destroy()
        deallocate(node)
     end do
  end do
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
  !![
  <objectDestructor name="accretionHalo_"/>
  !!]
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call nodeClassHierarchyFinalize         ()
  call parameters%reset  ()
  call parameters%destroy()
end program Test_Accretion_Halo_Metals

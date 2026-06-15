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

!!{RST
Contains a module which implements an N-body data operator which computes formation times of halos.
!!}

  use, intrinsic :: ISO_C_Binding      , only : c_size_t
  use            :: Cosmology_Functions, only : cosmologyFunctionsClass

  !![
  <nbodyOperator name="nbodyOperatorTimeOfFormation" docformat="rst">
   <description>
   An N-body data operator which computes formation times of halos.
   </description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorTimeOfFormation
     !!{RST
     An N-body data operator which computes formation times of halos.
     !!}
     private
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null()
     double precision                                   :: fractionFormation
   contains
     final     ::            timeOfFormationDestructor
     procedure :: operate => timeOfFormationOperate
  end type nbodyOperatorTimeOfFormation

  interface nbodyOperatorTimeOfFormation
     !!{RST
     Constructors for the "timeOfFormation" N-body operator class.
     !!}
     module procedure timeOfFormationConstructorParameters
     module procedure timeOfFormationConstructorInternal
  end interface nbodyOperatorTimeOfFormation

contains

  function timeOfFormationConstructorParameters(parameters) result (self)
    !!{RST
    Constructor for the "timeOfFormation" N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nbodyOperatorTimeOfFormation)                :: self
    type            (inputParameters             ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass     ), pointer       :: cosmologyFunctions_
    double precision                                              :: fractionFormation

    !![
    <inputParameter docformat="rst">
      <name>fractionFormation</name>
      <source>parameters</source>
      <description>
      The mass fraction used to define the formation epoch.
      </description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    self=nbodyOperatorTimeOfFormation(fractionFormation,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function timeOfFormationConstructorParameters

  function timeOfFormationConstructorInternal(fractionFormation,cosmologyFunctions_) result (self)
    !!{RST
    Internal constructor for the "timeOfFormation" N-body operator class.
    !!}
    implicit none
    type            (nbodyOperatorTimeOfFormation)                        :: self
    class           (cosmologyFunctionsClass     ), intent(in   ), target :: cosmologyFunctions_
    double precision                              , intent(in   )         :: fractionFormation
    !![
    <constructorAssign variables="fractionFormation, *cosmologyFunctions_"/>
    !!]

    return
  end function timeOfFormationConstructorInternal
  
  subroutine timeOfFormationDestructor(self)
    !!{RST
    Destructor for the "timeOfFormation" N-body operator class.
    !!}
    implicit none
    type(nbodyOperatorTimeOfFormation), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine timeOfFormationDestructor

  subroutine timeOfFormationOperate(self,simulations)
    !!{RST
    Compute formation times of halos.
    !!}
    use    :: Arrays_Search     , only : searchArray       , searchIndexed
    use    :: Display           , only : displayCounter    , displayCounterClear   , displayIndent, displayMessage, &
         &                               displayUnindent   , verbosityLevelStandard
    use    :: Error             , only : Error_Report
#ifdef USEMPI
    use    :: MPI_Utilities     , only : mpiSelf
#endif
    !$ use :: OMP_Lib           , only : OMP_Get_Thread_Num
    use    :: Sorting           , only : sortIndex
    implicit none
    class           (nbodyOperatorTimeOfFormation), intent(inout)                 :: self
    type            (nBodyData                   ), intent(inout), dimension(:  ) :: simulations
    double precision                              , pointer      , dimension(:  ) :: mass                    , expansionFactor, &
         &                                                                           timeFormation
    integer         (c_size_t                    ), pointer      , dimension(:  ) :: particleID              , descendantID
    integer         (c_size_t                    ), allocatable  , dimension(:  ) :: indexID
    class           (cosmologyFunctionsClass     ), pointer                       :: cosmologyFunctions_
    integer         (c_size_t                    )                                :: iSimulation             , i              , &
         &                                                                           j                       , l              , &
         &                                                                           m
    double precision                                                              :: expansionFactorFormation, massProgenitor
    
#ifdef USEMPI
    if (mpiSelf%isMaster()) then
#endif
       call displayIndent('compute formation times',verbosityLevelStandard)
#ifdef USEMPI
    end if
#endif
    ! Iterate over simulations.
    do iSimulation=1_c_size_t,size(simulations)
#ifdef USEMPI
       if (mpiSelf%isMaster()) then
#endif
       call displayMessage(var_str('simulation "')//simulations(iSimulation)%label//'"',verbosityLevelStandard)
#ifdef USEMPI
       end if
#endif
       ! Get the mass data.
       if (simulations(iSimulation)%propertiesReal   %exists('massVirial'     )) then
          mass            => simulations(iSimulation)%propertiesReal   %value('massVirial'     )
       else
          allocate(mass           (0))
          call Error_Report('halo virial masses are required, but are not available in the simulation'//{introspection:location})
       end if
       ! Get the expansion factor data.
       if (simulations(iSimulation)%propertiesReal   %exists('expansionFactor')) then
          expansionFactor => simulations(iSimulation)%propertiesReal   %value('expansionFactor')
       else
          allocate(expansionFactor(0))
          call Error_Report('expansion factors are required, but are not available in the simulation'//{introspection:location})
       end if
       ! Get the particle ID data.
       if (simulations(iSimulation)%propertiesInteger%exists('particleID'     )) then
          particleID      => simulations(iSimulation)%propertiesInteger%value('particleID'     )
       else
          allocate(particleID     (0))
          call Error_Report('particle IDs are required, but are not available in the simulation'     //{introspection:location})
       end if
       ! Get the descendant ID data.
       if (simulations(iSimulation)%propertiesInteger%exists('descendantID'   )) then
          descendantID    => simulations(iSimulation)%propertiesInteger%value('descendantID'   )
       else
          allocate(descendantID   (0))
          call Error_Report('descendant IDs are required, but are not available in the simulation'   //{introspection:location})
       end if
       ! Build an index into descendant IDs.
       indexID=sortIndex(descendantID)
       ! Allocate workspace.
       allocate(timeFormation(size(mass)))
       ! Compute formation times.
#ifdef USEMPI
       if (mpiSelf%isMaster()) then
#endif
          call displayCounter(0,.true.)
#ifdef USEMPI
       end if
#endif
       ! Visit all particles.
       !$omp parallel private(i,j,l,m,massProgenitor,expansionFactorFormation,cosmologyFunctions_)
       allocate(cosmologyFunctions_,mold=self%cosmologyFunctions_)
       !$omp critical(nBodyOperatorTimeFormationDeepCopy)
       !![
       <deepCopyReset variables="self%cosmologyFunctions_"/>
       <deepCopy source="self%cosmologyFunctions_" destination="cosmologyFunctions_"/>
       <deepCopyFinalize variables="cosmologyFunctions_"/>
       !!]
       !$omp end critical(nBodyOperatorTimeFormationDeepCopy)
       !$omp barrier
       !$omp do schedule(dynamic)
       do i=1_c_size_t,size(mass)
          ! Initialize the formation time to an impossible value.
          expansionFactorFormation=expansionFactor(i)
          ! Iterate over progenitor halos.
          j=i
          do while (j > 0_c_size_t)
             ! If this halo is above the formation threshold, and exists earlier than the current formation time, then update the formation time.
             if (mass(j) > self%fractionFormation*mass(i) .and. expansionFactor(j) < expansionFactorFormation) &
                  & expansionFactorFormation=expansionFactor(j)
             ! Find halos that descend into this halo.
             l=searchIndexed(descendantID,indexID,particleID(j))
             ! If no progenitors exist, exit the search.
             if (l < 1_c_size_t .or. l > size(descendantID)) exit
             ! Determine the most massive progenitor and its index.
             massProgenitor=-huge(0.0d0     )
             m             =-huge(0_c_size_t)
             do while (l > 0_c_size_t .and. descendantID(indexID(l)) == particleID(j))
                if (mass(indexID(l)) > massProgenitor) then
                   massProgenitor=mass(indexID(l))
                   m             =     indexID(l)
                end if
                l=l-1
             end do
             ! Move to the most massive progenitor.
             j=m             
          end do
          timeFormation(i)=cosmologyFunctions_%cosmicTime(expansionFactorFormation)
          ! Update progress.
          !$ if (OMP_Get_Thread_Num() == 0) then
#ifdef USEMPI
          if (mpiSelf%isMaster()) then
#endif
             call displayCounter(                        &
                  &              int(                    &
                  &                  +100.0d0            &
                  &                  *float(i         )  &
                  &                  /float(size(mass))  &
                  &                 )                  , &
                  &              .false.                 &
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
       call simulations(iSimulation)%propertiesReal%set('timeFormation',timeFormation)
       nullify(timeFormation)
    end do
#ifdef USEMPI
    if (mpiSelf%isMaster()) then
#endif
       call displayUnindent('done',verbosityLevelStandard)
#ifdef USEMPI
    end if
#endif
    return
  end subroutine timeOfFormationOperate


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

  use :: Cosmological_Density_Field    , only : cosmologicalMassVarianceClass
  use :: Cosmology_Functions           , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters          , only : cosmologyParametersClass
  use :: Excursion_Sets_Barriers       , only : excursionSetBarrierClass
  use :: Excursion_Sets_First_Crossings, only : excursionSetFirstCrossingClass
  use :: Halo_Mass_Functions           , only : haloMassFunctionClass
  use :: Power_Spectra                 , only : powerSpectrumClass

  !![
  <task name="taskExcursionSets">
   <description>A task which computes and outputs the halo mass function and related quantities.</description>
  </task>
  !!]
  type, extends(taskClass) :: taskExcursionSets
     !!{
     Implementation of a task which computes and outputs the halo mass function and related quantities.
     !!}
     private
     type            (varying_string                )          :: outputGroup
     integer                                                   :: massesPerDecade                      , timesPerDecade
     double precision                                          :: massMaximum                          , massMinimum    , &
          &                                                       timeMinimum                          , timeMaximum    , &
          &                                                       redshiftMinimum                      , redshiftMaximum
     logical                                                   :: nodeComponentsInitialized  =  .false.
     class           (cosmologyParametersClass      ), pointer :: cosmologyParameters_       => null()
     class           (cosmologyFunctionsClass       ), pointer :: cosmologyFunctions_        => null()
     class           (cosmologicalMassVarianceClass ), pointer :: cosmologicalMassVariance_  => null()
     class           (haloMassFunctionClass         ), pointer :: haloMassFunction_          => null()
     class           (excursionSetBarrierClass      ), pointer :: excursionSetBarrier_       => null()
     class           (excursionSetFirstCrossingClass), pointer :: excursionSetFirstCrossing_ => null()
     class           (powerSpectrumClass            ), pointer :: powerSpectrum_             => null()
   contains
     final     ::            excursionSetsDestructor
     procedure :: perform => excursionSetsPerform
  end type taskExcursionSets

  interface taskExcursionSets
     !!{
     Constructors for the \refClass{taskExcursionSets} task.
     !!}
     module procedure excursionSetsConstructorParameters
     module procedure excursionSetsConstructorInternal
  end interface taskExcursionSets

contains

  function excursionSetsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{taskExcursionSets} task class which takes a parameter set as input.
    !!}
    use :: Galacticus_Nodes, only : nodeClassHierarchyInitialize
    use :: Input_Parameters, only : inputParameter              , inputParameters
    use :: Node_Components , only : Node_Components_Initialize
    implicit none
    type            (taskExcursionSets             )                :: self
    type            (inputParameters               ), intent(inout) :: parameters
    class           (cosmologyParametersClass      ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass       ), pointer       :: cosmologyFunctions_
    class           (cosmologicalMassVarianceClass ), pointer       :: cosmologicalMassVariance_
    class           (haloMassFunctionClass         ), pointer       :: haloMassFunction_
    class           (excursionSetBarrierClass      ), pointer       :: excursionSetBarrier_
    class           (excursionSetFirstCrossingClass), pointer       :: excursionSetFirstCrossing_
    class           (powerSpectrumClass            ), pointer       :: powerSpectrum_
    type            (inputParameters               ), pointer       :: parametersRoot
    integer                                                         :: massesPerDecade           , timesPerDecade
    double precision                                                :: massMaximum               , massMinimum    , &
         &                                                             timeMinimum               , timeMaximum    , &
         &                                                             redshiftMinimum           , redshiftMaximum
    type            (varying_string                )                :: outputGroup

    ! Ensure the nodes objects are initialized.
    if (associated(parameters%parent)) then
       parametersRoot => parameters%parent
       do while (associated(parametersRoot%parent))
          parametersRoot => parametersRoot%parent
       end do
       call nodeClassHierarchyInitialize(parametersRoot)
       call Node_Components_Initialize  (parametersRoot)
    else
       parametersRoot => null()
       call nodeClassHierarchyInitialize(parameters    )
       call Node_Components_Initialize  (parameters    )
    end if
    self%nodeComponentsInitialized=.true.
    !![
    <inputParameter>
      <name>massMinimum</name>
      <defaultValue>1.0d10</defaultValue>
      <description>The minimum mass at which to tabulate excursion set solutions.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massMaximum</name>
      <defaultValue>1.0d15</defaultValue>
      <description>The maximum mass at which to tabulate excursion set solutions.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massesPerDecade</name>
      <defaultValue>10</defaultValue>
      <description>The number of points per decade of mass at which to tabulate excursion set solutions.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>redshiftMinimum</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The minimum redshift at which to tabulate excursion set solutions.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>redshiftMaximum</name>
      <defaultValue>10.0d0</defaultValue>
      <description>The maximum redshift at which to tabulate excursion set solutions.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>timesPerDecade</name>
      <defaultValue>10</defaultValue>
      <description>The number of points per decade of time at which to tabulate excursion set solutions.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>outputGroup</name>
      <defaultValue>var_str('excursionSets')</defaultValue>
      <description>The HDF5 output group within which to write excursion set solution data.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"       name="cosmologyParameters_"       source="parameters"/>
    <objectBuilder class="cosmologyFunctions"        name="cosmologyFunctions_"        source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance"  name="cosmologicalMassVariance_"  source="parameters"/>
    <objectBuilder class="haloMassFunction"          name="haloMassFunction_"          source="parameters"/>
    <objectBuilder class="excursionSetBarrier"       name="excursionSetBarrier_"       source="parameters"/>
    <objectBuilder class="excursionSetFirstCrossing" name="excursionSetFirstCrossing_" source="parameters"/>
    <objectBuilder class="powerSpectrum"             name="powerSpectrum_"             source="parameters"/>
    !!]
    ! Convert redshifts to times.
    timeMinimum=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftMaximum))
    timeMaximum=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftMinimum))
    self=taskExcursionSets(                            &
         &                 massMinimum               , &
         &                 massMaximum               , &
         &                 massesPerDecade           , &
         &                 timeMinimum               , &
         &                 timeMaximum               , &
         &                 timesPerDecade            , &
         &                 outputGroup               , &
         &                 cosmologyParameters_      , &
         &                 cosmologyFunctions_       , &
         &                 cosmologicalMassVariance_ , &
         &                 haloMassFunction_         , &
         &                 excursionSetBarrier_      , &
         &                 excursionSetFirstCrossing_, &
         &                 powerSpectrum_              &
         &                )
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"      />
    <objectDestructor name="cosmologyFunctions_"       />
    <objectDestructor name="cosmologicalMassVariance_" />
    <objectDestructor name="haloMassFunction_"         />
    <objectDestructor name="excursionSetBarrier_"      />
    <objectDestructor name="excursionSetFirstCrossing_"/>
    <objectDestructor name="powerSpectrum_"            />
    !!]
    return
  end function excursionSetsConstructorParameters

  function excursionSetsConstructorInternal(                            &
       &                                    massMinimum               , &
       &                                    massMaximum               , &
       &                                    massesPerDecade           , &
       &                                    timeMinimum               , &
       &                                    timeMaximum               , &
       &                                    timesPerDecade            , &
       &                                    outputGroup               , &
       &                                    cosmologyParameters_      , &
       &                                    cosmologyFunctions_       , &
       &                                    cosmologicalMassVariance_ , &
       &                                    haloMassFunction_         , &
       &                                    excursionSetBarrier_      , &
       &                                    excursionSetFirstCrossing_, &
       &                                    powerSpectrum_              &
       &                                   ) result(self)
    !!{
    Constructor for the \refClass{taskExcursionSets} task class which takes a parameter set as input.
    !!}
    implicit none
    type            (taskExcursionSets             )                        :: self
    class           (cosmologyParametersClass      ), intent(in   ), target :: cosmologyParameters_
    class           (cosmologyFunctionsClass       ), intent(in   ), target :: cosmologyFunctions_
    class           (cosmologicalMassVarianceClass ), intent(in   ), target :: cosmologicalMassVariance_
    class           (haloMassFunctionClass         ), intent(in   ), target :: haloMassFunction_
    class           (excursionSetBarrierClass      ), intent(in   ), target :: excursionSetBarrier_
    class           (excursionSetFirstCrossingClass), intent(in   ), target :: excursionSetFirstCrossing_
    class           (powerSpectrumClass            ), intent(in   ), target :: powerSpectrum_
    integer                                         , intent(in   )         :: massesPerDecade           , timesPerDecade
    double precision                                , intent(in   )         :: massMaximum               , massMinimum   , &
         &                                                                     timeMinimum               , timeMaximum
    type            (varying_string                ), intent(in   )         :: outputGroup
    !![
    <constructorAssign variables="massMinimum, massMaximum, massesPerDecade, timeMinimum, timeMaximum, timesPerDecade, outputGroup, *cosmologyParameters_, *cosmologyFunctions_, *cosmologicalMassVariance_, *haloMassFunction_, *excursionSetBarrier_, *excursionSetFirstCrossing_, *powerSpectrum_"/>
    !!]

    self%redshiftMinimum=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(timeMaximum))
    self%redshiftMaximum=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(timeMinimum))    
    return
  end function excursionSetsConstructorInternal

  subroutine excursionSetsDestructor(self)
    !!{
    Destructor for the \refClass{taskExcursionSets} task class.
    !!}
    use :: Node_Components, only : Node_Components_Uninitialize
    implicit none
    type(taskExcursionSets), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"      />
    <objectDestructor name="self%cosmologyFunctions_"       />
    <objectDestructor name="self%cosmologicalMassVariance_" />
    <objectDestructor name="self%haloMassFunction_"         />
    <objectDestructor name="self%excursionSetBarrier_"      />
    <objectDestructor name="self%excursionSetFirstCrossing_"/>
    <objectDestructor name="self%powerSpectrum_"            />
    !!]
    call self%outputGroup%destroy()
    if (self%nodeComponentsInitialized) call Node_Components_Uninitialize()
    return
  end subroutine excursionSetsDestructor

  subroutine excursionSetsPerform(self,status)
    !!{
    Compute and output the halo mass function.
    !!}
    use :: Display                 , only : displayIndent     , displayUnindent
    use :: Error                   , only : errorStatusSuccess
    use :: Output_HDF5             , only : outputFile
    use :: Galacticus_Nodes        , only : treeNode
    use :: IO_HDF5                 , only : hdf5Object
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Ranges        , only : Make_Range        , rangeTypeLogarithmic
    implicit none
    class           (taskExcursionSets), intent(inout), target           :: self
    integer                            , intent(  out), optional         :: status
    double precision                   , allocatable  , dimension(:    ) :: mass                    , time                    , &
         &                                                                  wavenumber
    double precision                   , allocatable  , dimension(:,:  ) :: barrier                 , firstCrossingProbability, &
         &                                                                  powerSpectrumValue      , variance                , &
         &                                                                  massFunctionDifferential
    double precision                   , allocatable  , dimension(:,:,:) :: firstCrossingRate
    integer                                                              :: iMass                   , jMass                   , &
         &                                                                  massCount               , timeCount               , &
         &                                                                  iTime
    type            (treeNode         )                                  :: node
    double precision                                                     :: varianceProgenitor
    type            (hdf5Object       )                                  :: outputGroup

    call displayIndent('Begin task: excursion sets')
#ifdef USEMPI
    ! Indicate that all MPI processes are coordinated in their work.
    call self%excursionSetFirstCrossing_%coordinatedMPI(.true.)
#endif
    ! Create grids of times and masses.
    timeCount=max(int(dble(self% timesPerDecade)*log10(self%timeMaximum/self%timeMinimum))+1,2)
    massCount=max(int(dble(self%massesPerDecade)*log10(self%massMaximum/self%massMinimum))+1,2)
    ! Allocate arrays.
    allocate(mass                    (massCount                    ))
    allocate(variance                (massCount          ,timeCount))
    allocate(wavenumber              (massCount                    ))
    allocate(powerSpectrumValue      (massCount          ,timeCount))
    allocate(time                    (                    timeCount))
    allocate(barrier                 (massCount          ,timeCount))
    allocate(firstCrossingProbability(massCount          ,timeCount))
    allocate(massFunctionDifferential(massCount          ,timeCount))
    allocate(firstCrossingRate       (massCount,massCount,timeCount))
    mass=Make_Range(self%massMinimum,self%massMaximum,massCount,rangeType=rangeTypeLogarithmic)
    time=Make_Range(self%timeMinimum,self%timeMaximum,timeCount,rangeType=rangeTypeLogarithmic)
    ! Set first crossing rates to unphysical values.
    firstCrossingRate=-1.0d0
    ! Iterate over times and masses.
    do iMass=1,massCount
       do iTime=1,timeCount
          if (iTime == 1) then
             wavenumber           (iMass      )=+(                                                                                    &
                  &                               +3.0d0                                                                              &
                  &                               *                                                           mass      (iMass      ) &
                  &                               /4.0d0                                                                              &
                  &                               /Pi                                                                                 &
                  &                               /self%cosmologyParameters_%densityCritical()                                        &
                  &                               /self%cosmologyParameters_%OmegaMatter    ()                                        &
                  &                              )**(-1.0d0/3.0d0)
          end if
          powerSpectrumValue      (iMass      ,iTime)=self%powerSpectrum_            %power       (wavenumber=wavenumber(iMass      )                                      ,time=time(iTime)                              )
          variance                (iMass      ,iTime)=self%cosmologicalMassVariance_ %rootVariance(mass      =mass      (iMass      )                                      ,time=time(iTime)                              )**2
          barrier                 (iMass      ,iTime)=self%excursionSetBarrier_      %barrier     (variance  =variance  (iMass,iTime)                                      ,time=time(iTime),node=node,rateCompute=.false.)
          firstCrossingProbability(iMass      ,iTime)=self%excursionSetFirstCrossing_%probability (variance  =variance  (iMass,iTime)                                      ,time=time(iTime),node=node                    )
          massFunctionDifferential(iMass      ,iTime)=self%haloMassFunction_         %differential(mass      =mass      (iMass      )                                      ,time=time(iTime),node=node                    )
          ! Compute halo branching rates.
          do jMass=1,iMass-1
             varianceProgenitor                      =self%cosmologicalMassVariance_ %rootVariance(mass      =mass      (jMass      )                                      ,time=time(iTime)                              )**2
             firstCrossingRate    (iMass,jMass,iTime)=self%excursionSetFirstCrossing_%rate        (variance  =variance  (iMass,iTime),varianceProgenitor=varianceProgenitor,time=time(iTime),node=node                    )
          end do
       end do
    end do
#ifdef USEMPI
    ! Indicate that all MPI processes are no longer coordinated in their work.
    call self%excursionSetFirstCrossing_%coordinatedMPI(.false.)
#endif
    ! Write results to the output file.
    outputGroup=outputFile%openGroup(char(self%outputGroup),'Group containing data relating to the excursion set problem.')
    call outputGroup%writeDataset(mass                    ,'mass'                    ,'The mass of the halo [M☉]'                       )
    call outputGroup%writeDataset(time                    ,'time'                    ,'The cosmic time [Gyr]'                           )
    call outputGroup%writeDataset(wavenumber              ,'wavenumber'              ,'The wavenumber associated with this mass [Mpc⁻¹]')
    call outputGroup%writeDataset(powerSpectrumValue      ,'powerSpectrum'           ,'The power spectrum at this wavenumber [Mpc³]'    )
    call outputGroup%writeDataset(variance                ,'variance'                ,'The variance on this mass scale'                 )
    call outputGroup%writeDataset(barrier                 ,'barrier'                 ,'The excursion set barrier for this variance'     )
    call outputGroup%writeDataset(firstCrossingProbability,'firstCrossingProbability','The first crossing probability'                  )
    call outputGroup%writeDataset(massFunctionDifferential,'massFunction'            ,'The halo mass function [Mpc⁻³ M☉⁻¹]'             )
    call outputGroup%writeDataset(firstCrossingRate       ,'firstCrossingRate'       ,'The first crossing rate [Gyr⁻¹]'                 )
    call outputGroup%close()
    ! Deallocate arrays.
    deallocate(mass                    )
    deallocate(variance                )
    deallocate(barrier                 )
    deallocate(firstCrossingProbability)
    deallocate(massFunctionDifferential)
    deallocate(wavenumber              )
    deallocate(powerSpectrumValue      )
    deallocate(firstCrossingRate       )
    if (present(status)) status=errorStatusSuccess
    call displayUnindent('Done task: excursion sets' )
    return
  end subroutine excursionSetsPerform

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

  use :: Cosmology_Functions        , only : cosmologyFunctionsClass
  use :: Cosmological_Velocity_Field, only : cosmologicalVelocityFieldClass
  use :: Dark_Matter_Halo_Scales    , only : darkMatterHaloScaleClass
  use :: Output_Times               , only : outputTimesClass

  !![
  <task name="taskVelocityField">
   <description>A task which computes and outputs the cosmological velocity field.</description>
  </task>
  !!]
  type, extends(taskClass) :: taskVelocityField
     !!{
     Implementation of a task which computes and outputs the cosmological velocity field.
     !!}
     private
     class           (cosmologyFunctionsClass       ), pointer :: cosmologyFunctions_        => null()
     class           (cosmologicalVelocityFieldClass), pointer :: cosmologicalVelocityField_ => null()
     class           (darkMatterHaloScaleClass      ), pointer :: darkMatterHaloScale_       => null()
     class           (outputTimesClass              ), pointer :: outputTimes_               => null()
     double precision                                          :: massMinimum                         , massMaximum
     integer                                                   :: pointsPerDecade
     type            (varying_string                )          :: outputGroup
   contains
     final     ::            velocityFieldDestructor
     procedure :: perform => velocityFieldPerform
  end type taskVelocityField

  interface taskVelocityField
     !!{
     Constructors for the \refClass{taskVelocityField} task.
     !!}
     module procedure velocityFieldConstructorParameters
     module procedure velocityFieldConstructorInternal
  end interface taskVelocityField

contains

  function velocityFieldConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{taskVelocityField} task class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (taskVelocityField               )                :: self
    type            (inputParameters                 ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass         ), pointer       :: cosmologyFunctions_
    class           (cosmologicalVelocityFieldClass  ), pointer       :: cosmologicalVelocityField_
    class           (darkMatterHaloScaleClass        ), pointer       :: darkMatterHaloScale_
    class           (outputTimesClass                ), pointer       :: outputTimes_
    double precision                                                  :: massMinimum               , massMaximum
    integer                                                           :: pointsPerDecade
    type            (varying_string                  )                :: outputGroup

    !![
    <inputParameter>
      <name>massMinimum</name>
      <defaultValue>1.0d10</defaultValue>
      <description>The minimum mass scale at which to tabulate the velocity field.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massMaximum</name>
      <defaultValue>1.0d+3</defaultValue>
      <description>The maximum mass scale at which to tabulate the velocity field.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>pointsPerDecade</name>
      <defaultValue>10</defaultValue>
      <description>The number of points per decade of mass at which to tabulate the velocity field.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>outputGroup</name>
      <defaultValue>var_str('.')</defaultValue>
      <description>The HDF5 output group within which to write velocity field data.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale"       name="darkMatterHaloScale_"       source="parameters"/>
    <objectBuilder class="cosmologyFunctions"        name="cosmologyFunctions_"        source="parameters"/>
    <objectBuilder class="cosmologicalVelocityField" name="cosmologicalVelocityField_" source="parameters"/>
    <objectBuilder class="outputTimes"               name="outputTimes_"               source="parameters"/>
    !!]
    self=taskVelocityField(massMinimum,massMaximum,pointsPerDecade,outputGroup,cosmologyFunctions_,cosmologicalVelocityField_,darkMatterHaloScale_,outputTimes_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"      />
    <objectDestructor name="cosmologyFunctions_"       />
    <objectDestructor name="cosmologicalVelocityField_"/>
    <objectDestructor name="outputTimes_"              />
    !!]
    return
  end function velocityFieldConstructorParameters

  function velocityFieldConstructorInternal(massMinimum,massMaximum,pointsPerDecade,outputGroup,cosmologyFunctions_,cosmologicalVelocityField_,darkMatterHaloScale_,outputTimes_) result(self)
    !!{
    Internal constructor for the \refClass{taskVelocityField} task class.
    !!}
    implicit none
    type            (taskVelocityField             )                        :: self
    class           (darkMatterHaloScaleClass      ), intent(in   ), target :: darkMatterHaloScale_
    class           (cosmologyFunctionsClass       ), intent(in   ), target :: cosmologyFunctions_
    class           (cosmologicalVelocityFieldClass), intent(in   ), target :: cosmologicalVelocityField_
    class           (outputTimesClass              ), intent(in   ), target :: outputTimes_
    double precision                                , intent(in   )         :: massMinimum               , massMaximum
    integer                                         , intent(in   )         :: pointsPerDecade
    type            (varying_string                ), intent(in   )         :: outputGroup
    !![
    <constructorAssign variables="massMinimum, massMaximum, pointsPerDecade, outputGroup, *cosmologyFunctions_, *cosmologicalVelocityField_, *darkMatterHaloScale_, *outputTimes_"/>
    !!]

    return
  end function velocityFieldConstructorInternal

  subroutine velocityFieldDestructor(self)
    !!{
    Destructor for the \refClass{taskVelocityField} task class.
    !!}
    implicit none
    type(taskVelocityField), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"      />
    <objectDestructor name="self%cosmologyFunctions_"       />
    <objectDestructor name="self%cosmologicalVelocityField_"/>
    <objectDestructor name="self%outputTimes_"              />
    !!]
    return
  end subroutine velocityFieldDestructor

  subroutine velocityFieldPerform(self,status)
    !!{
    Compute and output the halo mass function.
    !!}
    use            :: Display                         , only : displayIndent        , displayUnindent     , displayCounter, displayCounterClear, &
         &                                                     verbosityLevelWorking
    use            :: Calculations_Resets             , only : Calculations_Reset
    use            :: Error                           , only : errorStatusSuccess
    use            :: Output_HDF5                     , only : outputFile
    use            :: Galacticus_Nodes                , only : nodeComponentBasic   , treeNode
    use            :: IO_HDF5                         , only : hdf5Object
    use, intrinsic :: ISO_C_Binding                   , only : c_size_t
    use            :: Numerical_Constants_Astronomical, only : massSolar
    use            :: Numerical_Constants_Prefixes    , only : kilo
    use            :: Numerical_Ranges                , only : Make_Range           , rangeTypeLogarithmic
    use            :: String_Handling                 , only : operator(//)
    implicit none
    class           (taskVelocityField ), intent(inout), target           :: self
    integer                             , intent(  out), optional         :: status
    integer         (c_size_t          )                                  :: outputCount                     , massCount             , &
         &                                                                   iOutput                         , iMass                 , &
         &                                                                   jMass                           , countProgress         , &
         &                                                                   countTotal
    double precision                    , allocatable  , dimension(:    ) :: mass                            , epochTime             , &
         &                                                                   epochRedshift
    double precision                    , allocatable  , dimension(:,:  ) :: velocityDispersion1D
    double precision                    , allocatable  , dimension(:,:,:) :: velocityDispersion1DMergingHalos
    type            (treeNode          )               , pointer          :: node
    class           (nodeComponentBasic)               , pointer          :: basic
    type            (hdf5Object        )                                  :: outputsGroup                    , outputGroup           , &
         &                                                                   containerGroup                  , dataset
    type            (varying_string    )                                  :: groupName                       , description
    double precision                                                      :: radiusVirial                    , radiusVirialLagrangian

    call displayIndent('Begin task: velocity field')
    ! Get the requested output redshifts.
    outputCount      =self%outputTimes_%count()
    ! Compute number of tabulation points.
    massCount=int(log10(self%massMaximum/self%massMinimum)*dble(self%pointsPerDecade))+1
    ! Allocate arrays for velocity field.
    allocate(epochTime                       (outputCount))
    allocate(epochRedshift                   (outputCount))
    allocate(mass                            (massCount                      ))
    allocate(velocityDispersion1D            (massCount          ,outputCount))
    allocate(velocityDispersion1DMergingHalos(massCount,massCount,outputCount))
    ! Build a range of masses.
    mass(:)=Make_Range(self%massMinimum,self%massMaximum,int(massCount),rangeTypeLogarithmic)
    ! Construct a tree node.
    node  => treeNode      (                 )
    basic => node    %basic(autoCreate=.true.)
    ! Initialize progress counter.
    countTotal   =outputCount*(massCount*(massCount+1))/2
    countProgress=0    
    ! Iterate over outputs.
    do iOutput=1,outputCount
       epochTime    (iOutput)=                                                       self%outputTimes_%time(iOutput)
       epochRedshift(iOutput)=self%cosmologyFunctions_ %redshiftFromExpansionFactor(                                 &
            &                  self%cosmologyFunctions_%expansionFactor             (                                &
            &                                                                        self%outputTimes_%time(iOutput) &
            &                                                                       )                                &
            &                                                                      )
       call basic%timeSet(epochTime(iOutput))
       call basic%timeLastIsolatedSet(epochTime(iOutput))
       ! Iterate over all masses computing velocity field.
       do iMass=1,massCount
          ! Compute velocity dispersion.
          velocityDispersion1D(iMass,iOutput)=self%cosmologicalVelocityField_%velocityDispersion1D(time=self%outputTimes_%time(iOutput),mass=mass(iMass))
          ! Iterate over all masses.
          do jMass=1,iMass
             call displayCounter(int(100.0d0*dble(countProgress)/dble(countTotal)),isNew=countProgress==0,verbosity=verbosityLevelWorking)
             countProgress=countProgress+1          
              call basic%massSet(max(mass(iMass),mass(jMass)))
              call Calculations_Reset(node)
              radiusVirial                                         =+  self%darkMatterHaloScale_      %radiusVirial                    (           node                                             )
              radiusVirialLagrangian                               =+(                                                                                                                                &
                   &                                                  +self%darkMatterHaloScale_      %densityMean                     (           node                                             ) &
                   &                                                  /self%cosmologyFunctions_       %matterDensityEpochal            (           epochTime                               (iOutput)) &
                   &                                                 )**(1.0d0/3.0d0)                                                                                                                 &
                   &                                                *radiusVirial
              velocityDispersion1DMergingHalos(iMass,jMass,iOutput)=+  self%cosmologicalVelocityField_%velocityDispersion1DHaloPairwise(                                                              &
                   &                                                                                                                    time      =self%outputTimes_%time                  (iOutput), &
                   &                                                                                                                    mass1     =                  mass                  (iMass  ), &
                   &                                                                                                                    mass2     =                  mass                  (jMass  ), &
                   &                                                                                                                    separation=                  radiusVirialLagrangian           &
                   &                                                                                                                  )
              velocityDispersion1DMergingHalos(jMass,iMass,iOutput)=+velocityDispersion1DMergingHalos(iMass,jMass,iOutput)
           end do
       end do
    end do
    call displayCounterClear(verbosityLevelWorking)
    call node%destroy()
    deallocate(node)
    ! Open the group for output time information.
    if (self%outputGroup == ".") then
       outputsGroup  =outputFile    %openGroup(     'Outputs'        ,'Group containing datasets relating to output times.')
    else
       containerGroup=outputFile    %openGroup(char(self%outputGroup),'Group containing velocity field data.'              )
       outputsGroup  =containerGroup%openGroup(     'Outputs'        ,'Group containing datasets relating to output times.')
    end if
    ! Iterate over output times and output data.
    do iOutput=1,outputCount
       groupName  ='Output'
       description='Data for output number '
       groupName  =groupName  //iOutput
       description=description//iOutput
       outputGroup=outputsGroup%openGroup(char(groupName),char(description))
       call outputGroup   %writeAttribute(epochRedshift                   (    iOutput),'outputRedshift'                                                                                                   )
       call outputGroup   %writeAttribute(epochTime                       (    iOutput),'outputTime'                                                                                                       )
       call outputGroup   %writeDataset  (mass                                         ,'mass'                            ,'The mass.'                                             ,datasetReturned=dataset)
       call dataset       %writeAttribute(massSolar                                    ,'unitsInSI'                                                                                                        )
       call dataset       %close         (                                                                                                                                                                 )
       call outputGroup   %writeDataset  (velocityDispersion1D            (:  ,iOutput),'velocityDispersion1D'            ,'The 1-D velocity dispersion.'                          ,datasetReturned=dataset)
       call dataset       %writeAttribute(kilo                                         ,'unitsInSI'                                                                                                        )
       call dataset       %close         (                                                                                                                                                                 )
       call outputGroup   %writeDataset  (velocityDispersion1DMergingHalos(:,:,iOutput),'velocityDispersion1DMergingHalos','The 1-D velocity dispersion of pairs of merging halos.',datasetReturned=dataset)
       call dataset       %writeAttribute(kilo                                         ,'unitsInSI'                                                                                                        )
       call dataset       %close         (                                                                                                                                                                 )
        call outputGroup%close()
    end do
    call outputsGroup%close()
    if (self%outputGroup /= ".") call containerGroup%close()
    if (present(status)) status=errorStatusSuccess
    call displayUnindent('Done task: velocity field' )
    return
 end subroutine velocityFieldPerform

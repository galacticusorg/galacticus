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
  Implements an N-body data operator which computes mass functions.
  !!}
  
  use            :: Cosmology_Parameters, only : cosmologyParametersClass
  use, intrinsic :: ISO_C_Binding       , only : c_size_t

  !![
  <nbodyOperator name="nbodyOperatorSpinDistributionFunction">
   <description>An N-body data operator which computes mass functions.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorSpinDistributionFunction
     !!{
     An N-body data operator which computes mass functions.
     !!}
     private
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
     double precision                                    :: massMinimum                   , massMaximum  , &
          &                                                 spinMinimum                   , spinMaximum
     integer         (c_size_t                )          :: spinCountPerDecade
     type            (varying_string          )          :: simulationReference           , simulationURL, &
          &                                                 description
   contains
     final     ::            spinDistributionFunctionDestructor
     procedure :: operate => spinDistributionFunctionOperate
  end type nbodyOperatorSpinDistributionFunction

  interface nbodyOperatorSpinDistributionFunction
     !!{
     Constructors for the \refClass{nbodyOperatorSpinDistributionFunction} N-body operator class.
     !!}
     module procedure spinDistributionFunctionConstructorParameters
     module procedure spinDistributionFunctionConstructorInternal
  end interface nbodyOperatorSpinDistributionFunction

contains

  function spinDistributionFunctionConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{nbodyOperatorSpinDistributionFunction} N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nbodyOperatorSpinDistributionFunction)                :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    class           (cosmologyParametersClass             ), pointer       :: cosmologyParameters_
    double precision                                                       :: massMinimum         , massMaximum  , &
          &                                                                   spinMinimum         , spinMaximum
    integer         (c_size_t                             )                :: spinCountPerDecade
    type            (varying_string                       )                :: simulationReference , simulationURL, &
         &                                                                    description
    
    !![
    <inputParameter>
      <name>massMinimum</name>
      <source>parameters</source>
      <description>The minimum mass to consider.</description>
    </inputParameter>
    <inputParameter>
      <name>massMaximum</name>
      <source>parameters</source>
      <description>The maximum mass to consider.</description>
    </inputParameter>
    <inputParameter>
      <name>spinMinimum</name>
      <source>parameters</source>
      <description>The minimum spin to consider.</description>
    </inputParameter>
    <inputParameter>
      <name>spinMaximum</name>
      <source>parameters</source>
      <description>The maximum spin to consider.</description>
    </inputParameter>
    <inputParameter>
      <name>spinCountPerDecade</name>
      <source>parameters</source>
      <description>The number of bins per decade of spin.</description>
    </inputParameter>
    <inputParameter>
      <name>description</name>
      <source>parameters</source>
      <description>A description of this spin distribution function.</description>
    </inputParameter>
    <inputParameter>
      <name>simulationReference</name>
      <source>parameters</source>
      <description>A reference for the simulation.</description>
    </inputParameter>
    <inputParameter>
      <name>simulationURL</name>
      <source>parameters</source>
      <description>A URL for the simulation.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !!]
    self=nbodyOperatorSpinDistributionFunction(massMinimum,massMaximum,spinMinimum,spinMaximum,spinCountPerDecade,description,simulationReference,simulationURL,cosmologyParameters_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    !!]
    return
  end function spinDistributionFunctionConstructorParameters

  function spinDistributionFunctionConstructorInternal(massMinimum,massMaximum,spinMinimum,spinMaximum,spinCountPerDecade,description,simulationReference,simulationURL,cosmologyParameters_) result (self)
    !!{
    Internal constructor for the \refClass{nbodyOperatorSpinDistributionFunction} N-body operator class.
    !!}
    implicit none
    type            (nbodyOperatorSpinDistributionFunction)                        :: self
    double precision                                       , intent(in   )         :: massMinimum         , massMaximum, &
         &                                                                            spinMinimum         , spinMaximum
    integer         (c_size_t                             ), intent(in   )         :: spinCountPerDecade
    type            (varying_string                       ), intent(in   )         :: simulationReference , simulationURL, &
         &                                                                            description
    class           (cosmologyParametersClass             ), intent(in   ), target :: cosmologyParameters_
    !![
    <constructorAssign variables="massMinimum, massMaximum, spinMinimum, spinMaximum, spinCountPerDecade, description, simulationReference, simulationURL, *cosmologyParameters_"/>
    !!]
    
    return
  end function spinDistributionFunctionConstructorInternal
  
  subroutine spinDistributionFunctionDestructor(self)
    !!{
    Destructor for the \refClass{nbodyOperatorSpinDistributionFunction} N-body operator class.
    !!}
    implicit none
    type(nbodyOperatorSpinDistributionFunction), intent(inout) :: self
    
    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    !!]
    return
  end subroutine spinDistributionFunctionDestructor

  subroutine spinDistributionFunctionOperate(self,simulations)
    !!{
    Compute spin distribution function of particles.
    !!}
    use    :: Dates_and_Times   , only : Formatted_Date_and_Time
    use    :: Display           , only : displayCounter         , displayCounterClear   , displayIndent, displayMessage, &
         &                               displayUnindent        , verbosityLevelStandard
    use    :: Error             , only : Error_Report
    use    :: IO_HDF5           , only : hdf5Object
    use    :: HDF5_Access       , only : hdf5Access
    use    :: ISO_Varying_String, only : var_str
#ifdef USEMPI
    use    :: MPI_Utilities     , only : mpiSelf
#endif
    use    :: Numerical_Ranges  , only : Make_Range             , rangeTypeLogarithmic
    !$ use :: OMP_Lib           , only : OMP_Get_Thread_Num
    implicit none
    class           (nbodyOperatorSpinDistributionFunction), intent(inout)               :: self
    type            (nBodyData                            ), intent(inout), dimension(:) :: simulations
    double precision                                       , allocatable  , dimension(:) :: spinDistributionFunction, spinBin
    double precision                                       , pointer      , dimension(:) :: mass                    , spin
    integer         (c_size_t                             ), allocatable  , dimension(:) :: countBin
    integer         (c_size_t                             )                              :: iSimulation             , spinCount      , &
         &                                                                                  i                       , j
    integer                                                                              :: k
    double precision                                                                     :: binWidthInverse
    type            (hdf5Object                           )                              :: cosmologyGroup          , simulationGroup

#ifdef USEMPI
    if (mpiSelf%isMaster()) then
#endif
       call displayIndent('compute spin distribution function',verbosityLevelStandard)
#ifdef USEMPI
    end if
#endif
    ! Construct bins in mass.
    spinCount=int(log10(self%spinMaximum/self%spinMinimum)*dble(self%spinCountPerDecade),kind=c_size_t)
    allocate(spinBin                 (spinCount))
    allocate(spinDistributionFunction(spinCount))
    allocate(countBin                (spinCount))
    spinBin        =Make_Range(self%spinMinimum,self%spinMaximum,int(spinCount),rangeTypeLogarithmic,rangeBinned=.true.)
    binWidthInverse=1.0d0/log10(spinBin(2)/spinBin(1))
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
       if (simulations(iSimulation)%propertiesReal%exists('massVirial')) then
          mass => simulations(iSimulation)%propertiesReal%value('massVirial')
       else
          call Error_Report('halo virial masses are required, but are not available in the simulation'  //{introspection:location})
       end if
       ! Get the spin data.
       if (simulations(iSimulation)%propertiesReal%exists('spin')) then
          spin => simulations(iSimulation)%propertiesReal%value('spin')
       else
          call Error_Report('halo spin parameters are required, but are not available in the simulation'//{introspection:location})
       end if
       ! Accumulate counts.
#ifdef USEMPI
       if (mpiSelf%isMaster()) then
#endif
          call displayCounter(0,.true.)
#ifdef USEMPI
       end if
#endif
       countBin=0_c_size_t
       !$omp parallel do private(j) reduction(+:countBin) schedule(dynamic)
       do i=1_c_size_t,size(mass,kind=c_size_t)
#ifdef USEMPI
          ! If running under MPI with N processes, process only every Nth particle.
          if (mod(i,mpiSelf%count()) /= mpiSelf%rank()) cycle
#endif
          ! Skip particles not in the mass range.
          if     (                            &
               &   mass(i) < self%massMinimum &
               &  .or.                        &
               &   mass(i) > self%massMaximum &
               & ) cycle
          ! Accumulate particles into bins.
          j=int(log10(spin(i)/self%spinMinimum)*binWidthInverse)+1
          if (j >= 1 .and. j <= spinCount)  &
               & countBin(j)=+countBin  (j) &
               &             +1_c_size_t
          ! Update progress.
          !$ if (OMP_Get_Thread_Num() == 0) then
#ifdef USEMPI
          if (mpiSelf%isMaster()) then
#endif
             call displayCounter(                                      &
                  &              int(                                  &
                  &                  +100.0d0                          &
                  &                  *float(i                       )  &
                  &                  /float(size(mass,kind=c_size_t))  &
                  &                 )                                , &
                  &              .false.                               &
                  &             )
#ifdef USEMPI
          end if
#endif
          !$ end if
       end do
       !$omp end parallel do
#ifdef USEMPI
       if (mpiSelf%isMaster()) then
#endif
          call displayCounterClear()
#ifdef USEMPI
       end if
#endif
#ifdef USEMPI
       ! Reduce across MPI processes.
       countBin=mpiSelf%sum(countBin)
#endif
       ! Normalize the distribution function.
       spinDistributionFunction=dble(countBin)*binWidthInverse/log(10.0d0)/dble(sum(countBin))
#ifdef USEMPI
       if (mpiSelf%isMaster()) then
#endif
          call simulations(iSimulation)%analysis%writeDataset  (spinBin                  ,'spin'                    )
          call simulations(iSimulation)%analysis%writeDataset  (countBin                 ,'count'                   )
          call simulations(iSimulation)%analysis%writeDataset  (spinDistributionFunction ,'spinDistributionFunction')
          call simulations(iSimulation)%analysis%writeAttribute(self%description         ,"description"             )
          call simulations(iSimulation)%analysis%writeAttribute(Formatted_Date_and_Time(),"timestamp"               )
          cosmologyGroup=simulations(iSimulation)%analysis%openGroup('cosmology')
          call cosmologyGroup%writeAttribute(self%cosmologyParameters_%OmegaMatter    (),'OmegaMatter'    )
          call cosmologyGroup%writeAttribute(self%cosmologyParameters_%OmegaDarkEnergy(),'OmegaDarkEnergy')
          call cosmologyGroup%writeAttribute(self%cosmologyParameters_%HubbleConstant (),'HubbleConstant' )
          call cosmologyGroup%close         (                                                             )
          simulationGroup=simulations(iSimulation)%analysis%openGroup('simulation')
          call simulationGroup%writeAttribute(self%simulationReference,'reference')
          call simulationGroup%writeAttribute(self%simulationURL      ,'URL'      )
          do k=1,simulations(iSimulation)%attributesInteger%size()
             call simulationGroup%writeAttribute(simulations(iSimulation)%attributesInteger%value(k),char(simulations(iSimulation)%attributesInteger%key(k)))
          end do
          do k=1,simulations(iSimulation)%attributesReal   %size()
             call simulationGroup%writeAttribute(simulations(iSimulation)%attributesReal   %value(k),char(simulations(iSimulation)%attributesReal   %key(k)))
          end do
          do k=1,simulations(iSimulation)%attributesText   %size()
             call simulationGroup%writeAttribute(simulations(iSimulation)%attributesText   %value(k),char(simulations(iSimulation)%attributesText   %key(k)))
          end do
          call simulationGroup%close         (                                    )
          !$ call hdf5Access%unset()
#ifdef USEMPI
       end if
#endif
       nullify(mass)
       nullify(spin)
    end do
#ifdef USEMPI
    if (mpiSelf%isMaster()) then
#endif
       call displayUnindent('done',verbosityLevelStandard)
#ifdef USEMPI
    end if
#endif
    return
  end subroutine spinDistributionFunctionOperate


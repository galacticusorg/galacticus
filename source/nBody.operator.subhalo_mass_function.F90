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
Implements an N-body data operator which computes subhalo mass functions.
!!}

  use            :: Cosmology_Parameters, only : cosmologyParametersClass
  use, intrinsic :: ISO_C_Binding       , only : c_size_t

  !![
  <nbodyOperator name="nbodyOperatorSubhaloMassFunction">
   <description>An N-body data operator which computes subhalo mass functions.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorSubhaloMassFunction
     !!{
     An N-body data operator which computes subhalo mass functions.
     !!}
     private
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
     double precision                                    :: massRatioMinimum              , massRatioMaximum, &
          &                                                 massHost
     integer         (c_size_t                )          :: massCountPerDecade
     type            (varying_string          )          :: simulationReference           , simulationURL   , &
          &                                                 description
   contains
     final     ::            subhaloMassFunctionDestructor
     procedure :: operate => subhaloMassFunctionOperate
  end type nbodyOperatorSubhaloMassFunction

  interface nbodyOperatorSubhaloMassFunction
     !!{
     Constructors for the \refClass{nbodyOperatorSubhaloMassFunction} N-body operator class.
     !!}
     module procedure subhaloMassFunctionConstructorParameters
     module procedure subhaloMassFunctionConstructorInternal
  end interface nbodyOperatorSubhaloMassFunction

contains

  function subhaloMassFunctionConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{nbodyOperatorSubhaloMassFunction} N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nbodyOperatorSubhaloMassFunction)                :: self
    type            (inputParameters                 ), intent(inout) :: parameters
    class           (cosmologyParametersClass        ), pointer       :: cosmologyParameters_
    double precision                                                  :: massRatioMinimum    , massRatioMaximum, &
         &                                                               massHost
    integer         (c_size_t                        )                :: massCountPerDecade
    type            (varying_string                  )                :: simulationReference , simulationURL   , &
         &                                                               description

    !![
    <inputParameter>
      <name>massHost</name>
      <source>parameters</source>
      <description>The mass of the host halo.</description>
    </inputParameter>
    <inputParameter>
      <name>massRatioMinimum</name>
      <source>parameters</source>
      <description>The minimum mass ratio to consider.</description>
    </inputParameter>
    <inputParameter>
      <name>massRatioMaximum</name>
      <source>parameters</source>
      <description>The maximum mass ratio to consider.</description>
    </inputParameter>
    <inputParameter>
      <name>massCountPerDecade</name>
      <source>parameters</source>
      <description>The number of bins per decade of mass ratio.</description>
    </inputParameter>
    <inputParameter>
      <name>description</name>
      <source>parameters</source>
      <description>A description of this subhalo mass function.</description>
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
    self=nbodyOperatorSubhaloMassFunction(massHost,massRatioMinimum,massRatioMaximum,massCountPerDecade,description,simulationReference,simulationURL,cosmologyParameters_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    !!]
    return
  end function subhaloMassFunctionConstructorParameters

  function subhaloMassFunctionConstructorInternal(massHost,massRatioMinimum,massRatioMaximum,massCountPerDecade,description,simulationReference,simulationURL,cosmologyParameters_) result (self)
    !!{
    Internal constructor for the \refClass{nbodyOperatorSubhaloMassFunction} N-body operator class.
    !!}
    implicit none
    type            (nbodyOperatorSubhaloMassFunction)                        :: self
    double precision                                  , intent(in   )         :: massRatioMinimum    , massRatioMaximum, &
         &                                                                       massHost
    integer         (c_size_t                        ), intent(in   )         :: massCountPerDecade
    type            (varying_string                  ), intent(in   )         :: simulationReference , simulationURL, &
&                                                                                description
    class           (cosmologyParametersClass        ), intent(in   ), target :: cosmologyParameters_
    !![
    <constructorAssign variables="massHost, massRatioMinimum, massRatioMaximum, massCountPerDecade, description, simulationReference, simulationURL, *cosmologyParameters_"/>
    !!]

    return
  end function subhaloMassFunctionConstructorInternal
  
  subroutine subhaloMassFunctionDestructor(self)
    !!{
    Destructor for the \refClass{nbodyOperatorSubhaloMassFunction} N-body operator class.
    !!}
    implicit none
    type(nbodyOperatorSubhaloMassFunction), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    !!]
    return
  end subroutine subhaloMassFunctionDestructor

  subroutine subhaloMassFunctionOperate(self,simulations)
    !!{
    Compute mass functions of particles.
    !!}
    use    :: Dates_and_Times   , only : Formatted_Date_and_Time
    use    :: Display           , only : displayCounter         , displayCounterClear   , displayIndent, displayMessage, &
         &                               displayUnindent        , verbosityLevelStandard
    use    :: Error             , only : Error_Report
    use    :: HDF5_Access       , only : hdf5Access
    use    :: IO_HDF5           , only : hdf5Object
    use    :: ISO_Varying_String, only : var_str
#ifdef USEMPI
    use    :: MPI_Utilities     , only : mpiSelf
#endif
    use    :: Numerical_Ranges  , only : Make_Range             , rangeTypeLogarithmic
    !$ use :: OMP_Lib           , only : OMP_Get_Thread_Num
    implicit none
    class           (nbodyOperatorSubhaloMassFunction), intent(inout)               :: self
    type            (nBodyData                       ), intent(inout), dimension(:) :: simulations
    double precision                                  , allocatable  , dimension(:) :: massFunction     , massRatioBin
    double precision                                  , pointer      , dimension(:) :: mass
    integer         (c_size_t                        ), allocatable  , dimension(:) :: countBin
    integer         (c_size_t                        )                              :: iSimulation      , massCount         , &
         &                                                                             i                , j
    double precision                                                                :: binWidthInverse
    type            (hdf5Object                      )                              :: cosmologyGroup   , simulationGroup   , &
         &                                                                             massFunctionGroup

#ifdef USEMPI
    if (mpiSelf%isMaster()) then
#endif
       call displayIndent('compute subhalo mass function',verbosityLevelStandard)
#ifdef USEMPI
    end if
#endif
    ! Construct bins in mass.
    massCount=int(log10(self%massRatioMaximum/self%massRatioMinimum)*dble(self%massCountPerDecade),kind=c_size_t)
    allocate(massRatioBin     (massCount))
    allocate(massFunction(massCount))
    allocate(countBin    (massCount))
    massRatioBin   =Make_Range(self%massRatioMinimum,self%massRatioMaximum,int(massCount),rangeTypeLogarithmic,rangeBinned=.true.)
    binWidthInverse=1.0d0/log10(massRatioBin(2)/massRatioBin(1))
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
          call Error_Report('halo virial masses are required, but are not available in the simulation'//{introspection:location})
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
          ! Accumulate particles into bins.
          j=floor(log10(mass(i)/self%massHost/self%massRatioMinimum)*binWidthInverse)+1
          if (j >= 1 .and. j <= massCount)  &
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
       ! Compute mass function.
       massFunction=dble(countBin)*binWidthInverse/log(10.0d0)
#ifdef USEMPI
       if (mpiSelf%isMaster()) then
#endif
          massFunctionGroup=simulations(iSimulation)%analysis%openGroup('subhaloMassFunction')
          call massFunctionGroup%writeDataset  (massRatioBin             ,'massRatio'   )
          call massFunctionGroup%writeDataset  (countBin                 ,'count'       )
          call massFunctionGroup%writeDataset  (massFunction             ,'massFunction')
          call massFunctionGroup%writeAttribute(self%massHost            ,'massHost'    )
          call massFunctionGroup%writeAttribute(self%description         ,"description" )
          call massFunctionGroup%writeAttribute(Formatted_Date_and_Time(),"timestamp"   )
          call massFunctionGroup%close         (                                        )
          cosmologyGroup=simulations(iSimulation)%analysis%openGroup('cosmology')
          call cosmologyGroup%writeAttribute(self%cosmologyParameters_%OmegaMatter    (),'OmegaMatter'    )
          call cosmologyGroup%writeAttribute(self%cosmologyParameters_%OmegaDarkEnergy(),'OmegaDarkEnergy')
          call cosmologyGroup%writeAttribute(self%cosmologyParameters_%HubbleConstant (),'HubbleConstant' )
          call cosmologyGroup%close         (                                                             )
          simulationGroup=simulations(iSimulation)%analysis%openGroup('simulation')
          call simulationGroup%writeAttribute(self%simulationReference,'reference')
          call simulationGroup%writeAttribute(self%simulationURL      ,'URL'      )
          if (simulations(iSimulation)%attributesReal%exists('massParticle')) &
               & call simulationGroup%writeAttribute(simulations(iSimulation)%attributesReal%value('massParticle'),'massParticle')
          call simulationGroup%close         (                                    )
          !$ call hdf5Access%unset()
#ifdef USEMPI
       end if
#endif
       nullify(mass)
    end do
#ifdef USEMPI
    if (mpiSelf%isMaster()) then
#endif
       call displayUnindent('done',verbosityLevelStandard)
#ifdef USEMPI
    end if
#endif
    return
  end subroutine subhaloMassFunctionOperate


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
Implements an N-body data operator which computes subhalo radial distribution functions.
!!}

  use            :: Cosmology_Parameters, only : cosmologyParametersClass
  use, intrinsic :: ISO_C_Binding       , only : c_size_t

  !![
  <nbodyOperator name="nbodyOperatorSubhaloRadiusFunction">
   <description>An N-body data operator which computes subhalo radial distribution functions.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorSubhaloRadiusFunction
     !!{
     An N-body data operator which computes subhalo radial distribution functions.
     !!}
     private
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
     double precision                                    :: radiusRatioMinimum            , radiusRatioMaximum, &
          &                                                 radiusVirialHost              , massMinimum
     integer         (c_size_t                )          :: radiusCountPerDecade
     type            (varying_string          )          :: simulationReference           , simulationURL     , &
          &                                                 description
   contains
     final     ::            subhaloRadiusFunctionDestructor
     procedure :: operate => subhaloRadiusFunctionOperate
  end type nbodyOperatorSubhaloRadiusFunction

  interface nbodyOperatorSubhaloRadiusFunction
     !!{
     Constructors for the \refClass{nbodyOperatorSubhaloRadiusFunction} N-body operator class.
     !!}
     module procedure subhaloRadiusFunctionConstructorParameters
     module procedure subhaloRadiusFunctionConstructorInternal
  end interface nbodyOperatorSubhaloRadiusFunction

contains

  function subhaloRadiusFunctionConstructorParameters(parameters) result (self)
    !!{
    Constructor for the \refClass{nbodyOperatorSubhaloRadiusFunction} N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nbodyOperatorSubhaloRadiusFunction)                :: self
    type            (inputParameters                   ), intent(inout) :: parameters
    class           (cosmologyParametersClass          ), pointer       :: cosmologyParameters_
    double precision                                                    :: radiusRatioMinimum  , radiusRatioMaximum, &
         &                                                                 radiusVirialHost    , massMinimum
    integer         (c_size_t                          )                :: radiusCountPerDecade
    type            (varying_string                    )                :: simulationReference , simulationURL     , &
         &                                                                 description

    !![
    <inputParameter>
      <name>massMinimum</name>
      <source>parameters</source>
      <description>The minimum subhalo mass to include.</description>
    </inputParameter>
    <inputParameter>
      <name>radiusVirialHost</name>
      <source>parameters</source>
      <description>The virial radius of the host halo.</description>
    </inputParameter>
    <inputParameter>
      <name>radiusRatioMinimum</name>
      <source>parameters</source>
      <description>The minimum radius ratio to consider.</description>
    </inputParameter>
    <inputParameter>
      <name>radiusRatioMaximum</name>
      <source>parameters</source>
      <description>The maximum radius ratio to consider.</description>
    </inputParameter>
    <inputParameter>
      <name>radiusCountPerDecade</name>
      <source>parameters</source>
      <description>The number of bins per decade of radius ratio.</description>
    </inputParameter>
    <inputParameter>
      <name>description</name>
      <source>parameters</source>
      <description>A description of this subhalo radial distribution function.</description>
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
    self=nbodyOperatorSubhaloRadiusFunction(massMinimum,radiusVirialHost,radiusRatioMinimum,radiusRatioMaximum,radiusCountPerDecade,description,simulationReference,simulationURL,cosmologyParameters_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    !!]
    return
  end function subhaloRadiusFunctionConstructorParameters

  function subhaloRadiusFunctionConstructorInternal(massMinimum,radiusVirialHost,radiusRatioMinimum,radiusRatioMaximum,radiusCountPerDecade,description,simulationReference,simulationURL,cosmologyParameters_) result (self)
    !!{
    Internal constructor for the \refClass{nbodyOperatorSubhaloRadiusFunction} N-body operator class.
    !!}
    implicit none
    type            (nbodyOperatorSubhaloRadiusFunction)                        :: self
    double precision                                    , intent(in   )         :: radiusRatioMinimum  , radiusRatioMaximum, &
         &                                                                         radiusVirialHost    , massMinimum
    integer         (c_size_t                          ), intent(in   )         :: radiusCountPerDecade
    type            (varying_string                    ), intent(in   )         :: simulationReference , simulationURL     , &
&                                                                                  description
    class           (cosmologyParametersClass          ), intent(in   ), target :: cosmologyParameters_
    !![
    <constructorAssign variables="massMinimum, radiusVirialHost, radiusRatioMinimum, radiusRatioMaximum, radiusCountPerDecade, description, simulationReference, simulationURL, *cosmologyParameters_"/>
    !!]

    return
  end function subhaloRadiusFunctionConstructorInternal
  
  subroutine subhaloRadiusFunctionDestructor(self)
    !!{
    Destructor for the \refClass{nbodyOperatorSubhaloRadiusFunction} N-body operator class.
    !!}
    implicit none
    type(nbodyOperatorSubhaloRadiusFunction), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    !!]
    return
  end subroutine subhaloRadiusFunctionDestructor

  subroutine subhaloRadiusFunctionOperate(self,simulations)
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
    class           (nbodyOperatorSubhaloRadiusFunction), intent(inout)               :: self
    type            (nBodyData                         ), intent(inout), dimension(:) :: simulations
    double precision                                    , allocatable  , dimension(:) :: radialDistribution     , radiusRatioBin
    double precision                                    , pointer      , dimension(:) :: radius                 , mass
    integer         (c_size_t                          ), allocatable  , dimension(:) :: countBin
    integer         (c_size_t                          )                              :: iSimulation            , radiusCount    , &
         &                                                                               i                      , j
    double precision                                                                  :: binWidthInverse
    type            (hdf5Object                        )                              :: cosmologyGroup         , simulationGroup, &
         &                                                                               radialDistributionGroup

#ifdef USEMPI
    if (mpiSelf%isMaster()) then
#endif
       call displayIndent('compute subhalo radial distribution function',verbosityLevelStandard)
#ifdef USEMPI
    end if
#endif
    ! Construct bins in radius.
    radiusCount=int(log10(self%radiusRatioMaximum/self%radiusRatioMinimum)*dble(self%radiusCountPerDecade),kind=c_size_t)
    allocate(radiusRatioBin    (radiusCount))
    allocate(radialDistribution(radiusCount))
    allocate(countBin          (radiusCount))
    radiusRatioBin =Make_Range(self%radiusRatioMinimum,self%radiusRatioMaximum,int(radiusCount),rangeTypeLogarithmic,rangeBinned=.true.)
    binWidthInverse=1.0d0/log10(radiusRatioBin(2)/radiusRatioBin(1))
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
       if (simulations(iSimulation)%propertiesReal%exists('massVirial'       )) then
          mass   => simulations(iSimulation)%propertiesReal%value('massVirial'       )
       else
          call Error_Report('virial masses are required, but are not available in the simulation'                      //{introspection:location})
       end if
       ! Get the radius data.
       if (simulations(iSimulation)%propertiesReal%exists('distanceFromPoint')) then
          radius => simulations(iSimulation)%propertiesReal%value('distanceFromPoint')
       else
          call Error_Report('distances from the host halo center are required, but are not available in the simulation'//{introspection:location})
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
       do i=1_c_size_t,size(radius,kind=c_size_t)
#ifdef USEMPI
          ! If running under MPI with N processes, process only every Nth particle.
          if (mod(i,mpiSelf%count()) /= mpiSelf%rank        ()) cycle
#endif
          ! Skip halos below the mass limit.
          if (mass(i)                < self     %massMinimum  ) cycle
          ! Accumulate particles into bins.
          j=floor(log10(radius(i)/self%radiusVirialHost/self%radiusRatioMinimum)*binWidthInverse)+1
          if (j >= 1 .and. j <= radiusCount)  &
               & countBin(j)=+countBin  (j) &
               &             +1_c_size_t
          ! Update progress.
          !$ if (OMP_Get_Thread_Num() == 0) then
#ifdef USEMPI
          if (mpiSelf%isMaster()) then
#endif
             call displayCounter(                                        &
                  &              int(                                    &
                  &                  +100.0d0                            &
                  &                  *float(i                         )  &
                  &                  /float(size(radius,kind=c_size_t))  &
                  &                 )                                  , &
                  &              .false.                                 &
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
       ! Compute radial distribution function.
       radialDistribution=dble(countBin)*binWidthInverse/log(10.0d0)
#ifdef USEMPI
       if (mpiSelf%isMaster()) then
#endif
          radialDistributionGroup=simulations(iSimulation)%analysis%openGroup('subhaloRadiusFunction')
          call radialDistributionGroup%writeDataset  (radiusRatioBin           ,'radiusRatio'       )
          call radialDistributionGroup%writeDataset  (countBin                 ,'count'             )
          call radialDistributionGroup%writeDataset  (radialDistribution       ,'radialDistribution')
          call radialDistributionGroup%writeAttribute(self%massMinimum         ,'massMinimum'       )
          call radialDistributionGroup%writeAttribute(self%radiusVirialHost    ,'radiusVirialHost'  )
          call radialDistributionGroup%writeAttribute(self%description         ,"description"       )
          call radialDistributionGroup%writeAttribute(Formatted_Date_and_Time(),"timestamp"         )
          call radialDistributionGroup%close         (                                              )
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
       nullify(radius)
    end do
#ifdef USEMPI
    if (mpiSelf%isMaster()) then
#endif
       call displayUnindent('done',verbosityLevelStandard)
#ifdef USEMPI
    end if
#endif
    return
  end subroutine subhaloRadiusFunctionOperate


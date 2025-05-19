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
  <nbodyOperator name="nbodyOperatorConcentrationDistributionFunction">
   <description>An N-body data operator which computes mass functions.</description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorConcentrationDistributionFunction
     !!{
     An N-body data operator which computes mass functions.
     !!}
     private
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_        => null()
     double precision                                    :: massMinimum                          , massMaximum         , &
          &                                                 concentrationMinimum                 , concentrationMaximum
     integer         (c_size_t                )          :: concentrationCountPerDecade
     type            (varying_string          )          :: simulationReference                  , simulationURL       , &
          &                                                 description
   contains
     final     ::            concentrationDistributionFunctionDestructor
     procedure :: operate => concentrationDistributionFunctionOperate
  end type nbodyOperatorConcentrationDistributionFunction

  interface nbodyOperatorConcentrationDistributionFunction
     !!{
     Constructors for the {\normalfont \ttfamily concentrationDistributionFunction} N-body operator class.
     !!}
     module procedure concentrationDistributionFunctionConstructorParameters
     module procedure concentrationDistributionFunctionConstructorInternal
  end interface nbodyOperatorConcentrationDistributionFunction

contains

  function concentrationDistributionFunctionConstructorParameters(parameters) result (self)
    !!{
    Constructor for the {\normalfont \ttfamily concentrationDistributionFunction} N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nbodyOperatorConcentrationDistributionFunction)                :: self
    type            (inputParameters                               ), intent(inout) :: parameters
    class           (cosmologyParametersClass                      ), pointer       :: cosmologyParameters_
    double precision                                                                :: massMinimum                , massMaximum         , &
          &                                                                            concentrationMinimum       , concentrationMaximum
    integer         (c_size_t                                      )                :: concentrationCountPerDecade
    type            (varying_string                                )                :: simulationReference        , simulationURL       , &
         &                                                                             description

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
      <name>concentrationMinimum</name>
      <source>parameters</source>
      <description>The minimum concentration to consider.</description>
    </inputParameter>
    <inputParameter>
      <name>concentrationMaximum</name>
      <source>parameters</source>
      <description>The maximum concentration to consider.</description>
    </inputParameter>
    <inputParameter>
      <name>concentrationCountPerDecade</name>
      <source>parameters</source>
      <description>The number of bins per decade of concentration.</description>
    </inputParameter>
    <inputParameter>
      <name>description</name>
      <source>parameters</source>
      <description>A description of this concentration distribution function.</description>
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
    self=nbodyOperatorConcentrationDistributionFunction(massMinimum,massMaximum,concentrationMinimum,concentrationMaximum,concentrationCountPerDecade,description,simulationReference,simulationURL,cosmologyParameters_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    !!]
    return
  end function concentrationDistributionFunctionConstructorParameters

  function concentrationDistributionFunctionConstructorInternal(massMinimum,massMaximum,concentrationMinimum,concentrationMaximum,concentrationCountPerDecade,description,simulationReference,simulationURL,cosmologyParameters_) result (self)
    !!{
    Internal constructor for the {\normalfont \ttfamily concentrationDistributionFunction} N-body operator class.
    !!}
    implicit none
    type            (nbodyOperatorConcentrationDistributionFunction)                        :: self
    double precision                                                , intent(in   )         :: massMinimum                , massMaximum         , &
         &                                                                                     concentrationMinimum       , concentrationMaximum
    integer         (c_size_t                                      ), intent(in   )         :: concentrationCountPerDecade
    type            (varying_string                                ), intent(in   )         :: simulationReference        , simulationURL       , &
         &                                                                                     description
    class           (cosmologyParametersClass                      ), intent(in   ), target :: cosmologyParameters_
    !![
    <constructorAssign variables="massMinimum, massMaximum, concentrationMinimum, concentrationMaximum, concentrationCountPerDecade, description, simulationReference, simulationURL, *cosmologyParameters_"/>
    !!]
    
    return
  end function concentrationDistributionFunctionConstructorInternal
  
  subroutine concentrationDistributionFunctionDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily concentrationDistributionFunction} N-body operator class.
    !!}
    implicit none
    type(nbodyOperatorConcentrationDistributionFunction), intent(inout) :: self
    
    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    !!]
    return
  end subroutine concentrationDistributionFunctionDestructor

  subroutine concentrationDistributionFunctionOperate(self,simulations)
    !!{
    Compute concentration distribution function of particles.
    !!}
    use    :: Dates_and_Times   , only : Formatted_Date_and_Time
    use    :: Display           , only : displayCounter         , displayCounterClear   , displayIndent, displayMessage, &
         &                              displayUnindent        , verbosityLevelStandard
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
    class           (nbodyOperatorConcentrationDistributionFunction), intent(inout)               :: self
    type            (nBodyData                                     ), intent(inout), dimension(:) :: simulations
    double precision                                                , allocatable  , dimension(:) :: concentrationDistributionFunction, concentrationBin
    double precision                                                , pointer      , dimension(:) :: mass                             , radiusVirial      , &
         &                                                                                           radiusScale
    double precision                                                , allocatable  , dimension(:) :: concentration
    integer         (c_size_t                                      ), allocatable  , dimension(:) :: countBin
    integer         (c_size_t                                      )                              :: iSimulation                      , concentrationCount, &
         &                                                                                           i                                , j
    integer                                                                                       :: k
    double precision                                                                              :: binWidthInverse
    type            (hdf5Object                                    )                              :: cosmologyGroup                   , simulationGroup

#ifdef USEMPI
    if (mpiSelf%isMaster()) then
#endif
       call displayIndent('compute concentration distribution function',verbosityLevelStandard)
#ifdef USEMPI
    end if
#endif
    ! Construct bins in mass.
    concentrationCount=int(log10(self%concentrationMaximum/self%concentrationMinimum)*dble(self%concentrationCountPerDecade),kind=c_size_t)
    allocate(concentrationBin                 (concentrationCount))
    allocate(concentrationDistributionFunction(concentrationCount))
    allocate(countBin                         (concentrationCount))
    concentrationBin=Make_Range(self%concentrationMinimum,self%concentrationMaximum,int(concentrationCount),rangeTypeLogarithmic,rangeBinned=.true.)
    binWidthInverse =1.0d0/log10(concentrationBin(2)/concentrationBin(1))
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
       if (simulations(iSimulation)%propertiesReal%exists('massVirial'  )) then
          mass         => simulations(iSimulation)%propertiesReal%value('massVirial'  )
       else
          call Error_Report('halo virial masses are required, but are not available in the simulation'          //{introspection:location})
       end if
       ! Get the radii data.
       if (simulations(iSimulation)%propertiesReal%exists('radiusVirial')) then
          radiusVirial => simulations(iSimulation)%propertiesReal%value('radiusVirial')
       else
          allocate(radiusVirial(0))
          call Error_Report('halo virial radii parameters are required, but are not available in the simulation'//{introspection:location})
       end if
       if (simulations(iSimulation)%propertiesReal%exists('radiusScale' )) then
          radiusScale  => simulations(iSimulation)%propertiesReal%value('radiusScale' )
       else
          allocate(radiusScale (0))
          call Error_Report('halo scale radii parameters are required, but are not available in the simulation' //{introspection:location})
       end if
       ! Compute concentrations.
       allocate(concentration(size(mass,kind=c_size_t)))
       concentration=+radiusVirial &
            &        /radiusScale
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
          j=int(log10(concentration(i)/self%concentrationMinimum)*binWidthInverse)+1
          if (j >= 1 .and. j <= concentrationCount)  &
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
       concentrationDistributionFunction=dble(countBin)*binWidthInverse/log(10.0d0)/dble(sum(countBin))
#ifdef USEMPI
       if (mpiSelf%isMaster()) then
#endif
          call simulations(iSimulation)%analysis%writeDataset  (concentrationBin                 ,'concentration'                    )
          call simulations(iSimulation)%analysis%writeDataset  (countBin                         ,'count'                            )
          call simulations(iSimulation)%analysis%writeDataset  (concentrationDistributionFunction,'concentrationDistributionFunction')
          call simulations(iSimulation)%analysis%writeAttribute(self%description                 ,"description"                      )
          call simulations(iSimulation)%analysis%writeAttribute(Formatted_Date_and_Time()        ,"timestamp"                        )
          cosmologyGroup=simulations(iSimulation)%analysis%openGroup('cosmology')
          call cosmologyGroup%writeAttribute(self%cosmologyParameters_%OmegaMatter    (),'OmegaMatter'    )
          call cosmologyGroup%writeAttribute(self%cosmologyParameters_%OmegaDarkEnergy(),'OmegaDarkEnergy')
          call cosmologyGroup%writeAttribute(self%cosmologyParameters_%HubbleConstant (),'HubbleConstant' )
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
          !$ call hdf5Access%unset()
#ifdef USEMPI
       end if
#endif
       nullify   (mass         )
       nullify   (radiusVirial )
       nullify   (radiusScale  )
       deallocate(concentration)
    end do
#ifdef USEMPI
    if (mpiSelf%isMaster()) then
#endif
       call displayUnindent('done',verbosityLevelStandard)
#ifdef USEMPI
    end if
#endif
    return
  end subroutine concentrationDistributionFunctionOperate


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
  Implements an N-body data operator which computes mass functions.
  !!}
  
  use            :: Cosmology_Parameters, only : cosmologyParametersClass
  use, intrinsic :: ISO_C_Binding       , only : c_size_t

  !![
  <nbodyOperator name="nbodyOperatorSpinDistributionFunction" docformat="rst">
    <description>
    An N-body data operator which computes the halo spin distribution function by binning halos as a function of dimensionless spin parameter within a specified mass and spin range. Mass limits and binning are set by ``[massMinimum]``, ``[massMaximum]``, and ``[massCountPerDecade]``, spin limits and binning by ``[spinMinimum]``, ``[spinMaximum]``, and ``[spinCountPerDecade]``.
    </description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorSpinDistributionFunction
     !!{RST
     An N-body data operator which computes mass functions.
     !!}
     private
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
     double precision                                    :: massMinimum                   , massMaximum       , &
          &                                                 spinMinimum                   , spinMaximum
     integer         (c_size_t                )          :: spinCountPerDecade            , massCountPerDecade
     type            (varying_string          )          :: simulationReference           , simulationURL     , &
          &                                                 description
   contains
     final     ::            spinDistributionFunctionDestructor
     procedure :: operate => spinDistributionFunctionOperate
  end type nbodyOperatorSpinDistributionFunction

  interface nbodyOperatorSpinDistributionFunction
     !!{RST
     Constructors for the :galacticus-class:`nbodyOperatorSpinDistributionFunction` N-body operator class.
     !!}
     module procedure spinDistributionFunctionConstructorParameters
     module procedure spinDistributionFunctionConstructorInternal
  end interface nbodyOperatorSpinDistributionFunction

contains

  function spinDistributionFunctionConstructorParameters(parameters) result (self)
    !!{RST
    Constructor for the :galacticus-class:`nbodyOperatorSpinDistributionFunction` N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nbodyOperatorSpinDistributionFunction)                :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    class           (cosmologyParametersClass             ), pointer       :: cosmologyParameters_
    double precision                                                       :: massMinimum         , massMaximum       , &
          &                                                                   spinMinimum         , spinMaximum
    integer         (c_size_t                             )                :: spinCountPerDecade  , massCountPerDecade
    type            (varying_string                       )                :: simulationReference , simulationURL     , &
         &                                                                    description
    
    !![
    <inputParameter docformat="rst">
      <name>massMinimum</name>
      <source>parameters</source>
      <description>
      The minimum mass to consider.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>massMaximum</name>
      <source>parameters</source>
      <description>
      The maximum mass to consider.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>massCountPerDecade</name>
      <source>parameters</source>
      <description>
      The number of logarithmic bins per decade of mass used when constructing the spin distribution function.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>spinMinimum</name>
      <source>parameters</source>
      <description>
      The minimum dimensionless spin parameter below which halos are excluded from the spin distribution function histogram.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>spinMaximum</name>
      <source>parameters</source>
      <description>
      The maximum dimensionless spin parameter above which halos are excluded from the spin distribution function histogram.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>spinCountPerDecade</name>
      <source>parameters</source>
      <description>
      The number of logarithmic bins per decade of spin parameter used when constructing the spin distribution function.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>description</name>
      <source>parameters</source>
      <description>
      A human-readable description of this spin distribution function dataset, stored as metadata in the output file.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>simulationReference</name>
      <source>parameters</source>
      <description>
      A bibliographic reference for the N-body simulation from which this spin distribution is derived, stored as output metadata.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>simulationURL</name>
      <source>parameters</source>
      <description>
      A URL pointing to the publicly accessible dataset or documentation for the N-body simulation, stored as output metadata.
      </description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !!]
    self=nbodyOperatorSpinDistributionFunction(massMinimum,massMaximum,massCountPerDecade,spinMinimum,spinMaximum,spinCountPerDecade,description,simulationReference,simulationURL,cosmologyParameters_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    !!]
    return
  end function spinDistributionFunctionConstructorParameters

  function spinDistributionFunctionConstructorInternal(massMinimum,massMaximum,massCountPerDecade,spinMinimum,spinMaximum,spinCountPerDecade,description,simulationReference,simulationURL,cosmologyParameters_) result (self)
    !!{RST
    Internal constructor for the :galacticus-class:`nbodyOperatorSpinDistributionFunction` N-body operator class.
    !!}
    implicit none
    type            (nbodyOperatorSpinDistributionFunction)                        :: self
    double precision                                       , intent(in   )         :: massMinimum         , massMaximum       , &
         &                                                                            spinMinimum         , spinMaximum
    integer         (c_size_t                             ), intent(in   )         :: spinCountPerDecade  , massCountPerDecade
    type            (varying_string                       ), intent(in   )         :: simulationReference , simulationURL     , &
         &                                                                            description
    class           (cosmologyParametersClass             ), intent(in   ), target :: cosmologyParameters_
    !![
    <constructorAssign variables="massMinimum, massMaximum, massCountPerDecade, spinMinimum, spinMaximum, spinCountPerDecade, description, simulationReference, simulationURL, *cosmologyParameters_"/>
    !!]
    
    return
  end function spinDistributionFunctionConstructorInternal
  
  subroutine spinDistributionFunctionDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`nbodyOperatorSpinDistributionFunction` N-body operator class.
    !!}
    implicit none
    type(nbodyOperatorSpinDistributionFunction), intent(inout) :: self
    
    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    !!]
    return
  end subroutine spinDistributionFunctionDestructor

  subroutine spinDistributionFunctionOperate(self,simulations)
    !!{RST
    Compute spin distribution function of particles.
    !!}
    use    :: Dates_and_Times   , only : Formatted_Date_and_Time
    use    :: Display           , only : displayCounter         , displayCounterClear   , displayIndent, displayMessage, &
         &                               displayUnindent        , verbosityLevelStandard
    use    :: Error             , only : Error_Report
    use    :: IO_HDF5           , only : hdf5Group
    use    :: HDF5_Access       , only : hdf5Access
    use    :: ISO_Varying_String, only : var_str
#ifdef USEMPI
    use    :: MPI_Utilities     , only : mpiSelf
#endif
    use    :: Numerical_Ranges  , only : Make_Range             , rangeTypeLogarithmic
    !$ use :: OMP_Lib           , only : OMP_Get_Thread_Num
    implicit none
    class           (nbodyOperatorSpinDistributionFunction), intent(inout)                 :: self
    type            (nBodyData                            ), intent(inout), dimension(:  ) :: simulations
    double precision                                       , allocatable  , dimension(:  ) :: massBin, spinBin
    double precision                                       , allocatable  , dimension(:,:) :: spinDistributionFunction
    double precision                                       , pointer      , dimension(:  ) :: mass                    , spin
    integer         (c_size_t                             ), allocatable  , dimension(:,:) :: countBin
    integer         (c_size_t                             )                                :: iSimulation             , spinCount      , &
         &                                                                                    i                       , j              , &
         &                                                                                    k                       , massCount
    integer                                                                                :: m
    double precision                                                                       :: binWidthInverseSpin     , binWidthInverseMass, &
         &                                                                                    countTotal
    type            (hdf5Group                            )                                :: cosmologyGroup          , simulationGroup

#ifdef USEMPI
    if (mpiSelf%isMaster()) then
#endif
       call displayIndent('compute spin distribution function',verbosityLevelStandard)
#ifdef USEMPI
    end if
#endif
    ! Construct bins in mass.
    massCount=int(log10(self%massMaximum/self%massMinimum)*dble(self%massCountPerDecade),kind=c_size_t)
    spinCount=int(log10(self%spinMaximum/self%spinMinimum)*dble(self%spinCountPerDecade),kind=c_size_t)
    allocate(massBin                 (massCount          ))
    allocate(spinBin                 (          spinCount))
    allocate(spinDistributionFunction(massCount,spinCount))
    allocate(countBin                (massCount,spinCount))
    massBin            =Make_Range(self%massMinimum,self%massMaximum,int(massCount),rangeTypeLogarithmic,rangeBinned=.true.)
    spinBin            =Make_Range(self%spinMinimum,self%spinMaximum,int(spinCount),rangeTypeLogarithmic,rangeBinned=.true.)
    binWidthInverseMass=1.0d0/log10(massBin(2)/massBin(1))
    binWidthInverseSpin=1.0d0/log10(spinBin(2)/spinBin(1))
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
          j=floor(log10(mass(i)/self%massMinimum)*binWidthInverseMass)+1
          k=floor(log10(spin(i)/self%spinMinimum)*binWidthInverseSpin)+1
          if     (                              &
               &   j >= 1 .and. j <= massCount  &
               &  .and.                         &
               &   k >= 1 .and. k <= spinCount) &
               & countBin(j,k)=+countBin  (j,k) &
               &               +1_c_size_t
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
       do i=1_c_size_t,massCount
          countTotal                      =dble(sum(countBin(i,:)))
          if (countTotal > 0.0d0)  then
             spinDistributionFunction(i,:)=dble(    countBin(i,:) )*binWidthInverseSpin/log(10.0d0)/countTotal
          else
             spinDistributionFunction(i,:)=0.0d0
          end if
       end do
#ifdef USEMPI
       if (mpiSelf%isMaster()) then
#endif
          !$ call hdf5Access%set()
          call simulations(iSimulation)%analysis%writeDataset  (massBin                  ,'mass'                    )
          call simulations(iSimulation)%analysis%writeDataset  (spinBin                  ,'spin'                    )
          call simulations(iSimulation)%analysis%writeDataset  (countBin                 ,'count'                   )
          call simulations(iSimulation)%analysis%writeDataset  (spinDistributionFunction ,'spinDistributionFunction')
          call simulations(iSimulation)%analysis%writeAttribute(self%description         ,"description"             )
          call simulations(iSimulation)%analysis%writeAttribute(Formatted_Date_and_Time(),"timestamp"               )
          cosmologyGroup=simulations(iSimulation)%analysis%openGroup('cosmology')
          call cosmologyGroup%writeAttribute(self%cosmologyParameters_%OmegaMatter    (),'OmegaMatter'    )
          call cosmologyGroup%writeAttribute(self%cosmologyParameters_%OmegaDarkEnergy(),'OmegaDarkEnergy')
          call cosmologyGroup%writeAttribute(self%cosmologyParameters_%HubbleConstant (),'HubbleConstant' )
          simulationGroup=simulations(iSimulation)%analysis%openGroup('simulation')
          call simulationGroup%writeAttribute(self%simulationReference,'reference')
          call simulationGroup%writeAttribute(self%simulationURL      ,'URL'      )
          do m=1,simulations(iSimulation)%attributesInteger%size()
             call simulationGroup%writeAttribute(simulations(iSimulation)%attributesInteger%value(m),char(simulations(iSimulation)%attributesInteger%key(m)))
          end do
          do m=1,simulations(iSimulation)%attributesReal   %size()
             call simulationGroup%writeAttribute(simulations(iSimulation)%attributesReal   %value(m),char(simulations(iSimulation)%attributesReal   %key(m)))
          end do
          do m=1,simulations(iSimulation)%attributesText   %size()
             call simulationGroup%writeAttribute(simulations(iSimulation)%attributesText   %value(m),char(simulations(iSimulation)%attributesText   %key(m)))
          end do
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


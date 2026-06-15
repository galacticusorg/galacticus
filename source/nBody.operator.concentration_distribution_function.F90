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
  <nbodyOperator name="nbodyOperatorConcentrationDistributionFunction" docformat="rst">
    <description>
    An N-body data operator which computes the halo concentration distribution function by binning halos as a function of concentration within a specified mass and concentration range. Mass limits and binning are set by ``[massMinimum]`` , ``[massMaximum]``, and ``[massCountPerDecade]``, concentration limits and binning by ``[concentrationMinimum]``, ``[concentrationMaximum]``, and ``[concentrationCountPerDecade]``.
    </description>
  </nbodyOperator>
  !!]
  type, extends(nbodyOperatorClass) :: nbodyOperatorConcentrationDistributionFunction
     !!{RST
     An N-body data operator which computes mass functions.
     !!}
     private
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_        => null()
     double precision                                    :: massMinimum                          , massMaximum         , &
          &                                                 concentrationMinimum                 , concentrationMaximum
     integer         (c_size_t                )          :: concentrationCountPerDecade          , massCountPerDecade
     type            (varying_string          )          :: simulationReference                  , simulationURL       , &
          &                                                 description
   contains
     final     ::            concentrationDistributionFunctionDestructor
     procedure :: operate => concentrationDistributionFunctionOperate
  end type nbodyOperatorConcentrationDistributionFunction

  interface nbodyOperatorConcentrationDistributionFunction
     !!{RST
     Constructors for the :galacticus-class:`nbodyOperatorConcentrationDistributionFunction` N-body operator class.
     !!}
     module procedure concentrationDistributionFunctionConstructorParameters
     module procedure concentrationDistributionFunctionConstructorInternal
  end interface nbodyOperatorConcentrationDistributionFunction

contains

  function concentrationDistributionFunctionConstructorParameters(parameters) result (self)
    !!{RST
    Constructor for the :galacticus-class:`nbodyOperatorConcentrationDistributionFunction` N-body operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nbodyOperatorConcentrationDistributionFunction)                :: self
    type            (inputParameters                               ), intent(inout) :: parameters
    class           (cosmologyParametersClass                      ), pointer       :: cosmologyParameters_
    double precision                                                                :: massMinimum                , massMaximum         , &
          &                                                                            concentrationMinimum       , concentrationMaximum
    integer         (c_size_t                                      )                :: concentrationCountPerDecade, massCountPerDecade
    type            (varying_string                                )                :: simulationReference        , simulationURL       , &
         &                                                                             description

    !![
    <inputParameter docformat="rst">
      <name>massMinimum</name>
      <source>parameters</source>
      <description>
      The minimum halo mass (in :math:`\mathrm{M}_\odot`) below which halos are excluded from the concentration distribution function.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>massMaximum</name>
      <source>parameters</source>
      <description>
      The maximum halo mass (in :math:`\mathrm{M}_\odot`) above which halos are excluded from the concentration distribution function.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>massCountPerDecade</name>
      <source>parameters</source>
      <description>
      The number of logarithmic bins per decade of mass used when constructing the concentration distribution function.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>concentrationMinimum</name>
      <source>parameters</source>
      <description>
      The minimum halo concentration parameter below which halos are excluded from the distribution function histogram.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>concentrationMaximum</name>
      <source>parameters</source>
      <description>
      The maximum halo concentration parameter above which halos are excluded from the distribution function histogram.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>concentrationCountPerDecade</name>
      <source>parameters</source>
      <description>
      The number of logarithmic bins per decade of concentration parameter used when constructing the concentration distribution function.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>description</name>
      <source>parameters</source>
      <description>
      A human-readable description of this concentration distribution function dataset, stored as metadata in the output file.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>simulationReference</name>
      <source>parameters</source>
      <description>
      A bibliographic reference for the N-body simulation from which this concentration distribution is derived, stored as output metadata.
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
    self=nbodyOperatorConcentrationDistributionFunction(massMinimum,massMaximum,massCountPerDecade,concentrationMinimum,concentrationMaximum,concentrationCountPerDecade,description,simulationReference,simulationURL,cosmologyParameters_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    !!]
    return
  end function concentrationDistributionFunctionConstructorParameters

  function concentrationDistributionFunctionConstructorInternal(massMinimum,massMaximum,massCountPerDecade,concentrationMinimum,concentrationMaximum,concentrationCountPerDecade,description,simulationReference,simulationURL,cosmologyParameters_) result (self)
    !!{RST
    Internal constructor for the :galacticus-class:`nbodyOperatorConcentrationDistributionFunction` N-body operator class.
    !!}
    implicit none
    type            (nbodyOperatorConcentrationDistributionFunction)                        :: self
    double precision                                                , intent(in   )         :: massMinimum                , massMaximum         , &
         &                                                                                     concentrationMinimum       , concentrationMaximum
    integer         (c_size_t                                      ), intent(in   )         :: concentrationCountPerDecade, massCountPerDecade
    type            (varying_string                                ), intent(in   )         :: simulationReference        , simulationURL       , &
         &                                                                                     description
    class           (cosmologyParametersClass                      ), intent(in   ), target :: cosmologyParameters_
    !![
    <constructorAssign variables="massMinimum, massMaximum, massCountPerDecade, concentrationMinimum, concentrationMaximum, concentrationCountPerDecade, description, simulationReference, simulationURL, *cosmologyParameters_"/>
    !!]
    
    return
  end function concentrationDistributionFunctionConstructorInternal
  
  subroutine concentrationDistributionFunctionDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`nbodyOperatorConcentrationDistributionFunction` N-body operator class.
    !!}
    implicit none
    type(nbodyOperatorConcentrationDistributionFunction), intent(inout) :: self
    
    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    !!]
    return
  end subroutine concentrationDistributionFunctionDestructor

  subroutine concentrationDistributionFunctionOperate(self,simulations)
    !!{RST
    Compute concentration distribution function of particles.
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
    class           (nbodyOperatorConcentrationDistributionFunction), intent(inout)                 :: self
    type            (nBodyData                                     ), intent(inout), dimension(:  ) :: simulations
    double precision                                                , allocatable  , dimension(:  ) :: concentrationBin                 , concentration      , &
         &                                                                                             massBin
    double precision                                                , allocatable  , dimension(:,:) :: concentrationDistributionFunction
    double precision                                                , pointer      , dimension(:  ) :: mass                             , radiusVirial       , &
         &                                                                                             radiusScale
    integer         (c_size_t                                      ), allocatable  , dimension(:,:) :: countBin
    integer         (c_size_t                                      )                                :: iSimulation                      , concentrationCount , &
         &                                                                                             i                                , j                  , &
         &                                                                                             k                                , massCount
    integer                                                                                         :: m
    double precision                                                                                :: binWidthInverseConcentration     , binWidthInverseMass, &
         &                                                                                             countTotal
    type            (hdf5Object                                    )                                :: cosmologyGroup                   , simulationGroup

#ifdef USEMPI
    if (mpiSelf%isMaster()) then
#endif
       call displayIndent('compute concentration distribution function',verbosityLevelStandard)
#ifdef USEMPI
    end if
#endif
    ! Construct bins in mass and concentration.
    massCount         =int(log10(self%massMaximum         /self%massMinimum         )*dble(self%massCountPerDecade         ),kind=c_size_t)
    concentrationCount=int(log10(self%concentrationMaximum/self%concentrationMinimum)*dble(self%concentrationCountPerDecade),kind=c_size_t)
    allocate(massBin                          (massCount                   ))
    allocate(concentrationBin                 (          concentrationCount))
    allocate(concentrationDistributionFunction(massCount,concentrationCount))
    allocate(countBin                         (massCount,concentrationCount))
    massBin                     =Make_Range(self%massMinimum         ,self%massMaximum         ,int(massCount         ),rangeTypeLogarithmic,rangeBinned=.true.)
    concentrationBin            =Make_Range(self%concentrationMinimum,self%concentrationMaximum,int(concentrationCount),rangeTypeLogarithmic,rangeBinned=.true.)
    binWidthInverseMass         =1.0d0/log10(massBin         (2)/massBin         (1))
    binWidthInverseConcentration=1.0d0/log10(concentrationBin(2)/concentrationBin(1))
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
          ! Skip halos not in the mass range.
          if     (                            &
               &   mass(i) < self%massMinimum &
               &  .or.                        &
               &   mass(i) > self%massMaximum &
               & ) cycle
          ! Accumulate particles into bins.
          j=floor(log10(mass         (i)/self%massMinimum         )*binWidthInverseMass         )+1
          k=floor(log10(concentration(i)/self%concentrationMinimum)*binWidthInverseConcentration)+1
          if     (                                      &
               &   j >= 1 .and. j <= massCount          &
               &  .and.                                 &
               &   k >= 1 .and. k <= concentrationCount &
               & )                                      &
               & countBin(j,k)=+countBin  (j,k)         &
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
          countTotal                               =dble(sum(countBin(i,:)))
          if (countTotal > 0.0d0) then
             concentrationDistributionFunction(i,:)=dble(    countBin(i,:) )*binWidthInverseConcentration/log(10.0d0)/countTotal
          else
             concentrationDistributionFunction(i,:)=0.0d0
          end if
       end do
#ifdef USEMPI
       if (mpiSelf%isMaster()) then
#endif
          !$ call hdf5Access%set()
          call simulations(iSimulation)%analysis%writeDataset  (massBin                          ,'mass'                             )
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


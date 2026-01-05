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

  !+    Contributions to this file made by: Andrew Benson, Xiaolong Du.

  !!{
  Implements an output analysis class that computes subhalo $V_\mathrm{max}$ vs. $M_\mathrm{bound}$ tidal tracks.
  !!}

  use    :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass
  !$ use :: OMP_Lib                 , only : omp_lock_kind

  !![
  <outputAnalysis name="outputAnalysisTidalTracksVelocityMaximum">
    <description>An output analysis class that computes subhalo mean maximum velocity as a function of mass.</description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisClass) :: outputAnalysisTidalTracksVelocityMaximum
     !!{
     An output analysis class that computes subhalo mean maximum velocity as a function of mass.
     !!}
     private
     class           (darkMatterProfileDMOClass), pointer                   :: darkMatterProfileDMO_         => null(), darkMatterProfileDMOUnheated => null()
     !$ integer      (omp_lock_kind            )                            :: accumulateLock
     integer                                                                :: countPointsOnTrack
     double precision                                                       :: mu                                     , eta
     double precision                                                       :: logLikelihood_
     double precision                           , dimension(:), allocatable :: fractionMassBound                      , fractionVelocityMaximum               , &
          &                                                                    fractionVelocityMaximumTarget
   contains
     final     ::                  tidalTracksVelocityMaximumDestructor
     procedure :: analyze       => tidalTracksVelocityMaximumAnalyze
     procedure :: finalize      => tidalTracksVelocityMaximumFinalize
     procedure :: reduce        => tidalTracksVelocityMaximumReduce
     procedure :: logLikelihood => tidalTracksVelocityMaximumLogLikelihood
  end type outputAnalysisTidalTracksVelocityMaximum

  interface outputAnalysisTidalTracksVelocityMaximum
     !!{
     Constructors for the \refClass{outputAnalysisTidalTracksVelocityMaximum} output analysis class.
     !!}
     module procedure tidalTracksVelocityMaximumConstructorParameters
     module procedure tidalTracksVelocityMaximumConstructorInternal
  end interface outputAnalysisTidalTracksVelocityMaximum

contains

  function tidalTracksVelocityMaximumConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisTidalTracksVelocityMaximum} output analysis class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (outputAnalysisTidalTracksVelocityMaximum)                :: self
    type            (inputParameters                         ), intent(inout) :: parameters
    class           (darkMatterProfileDMOClass               ), pointer       :: darkMatterProfileDMO_, darkMatterProfileDMOUnheated
    double precision                                                          :: mu                   , eta
    
    !![
    <inputParameter>
      <name>mu</name>
      <source>parameters</source>
      <defaultValue>0.4d0</defaultValue>
      <description>The parameter $\mu$ in the \cite{penarrubia_impact_2010} tidal track fitting function.</description>
    </inputParameter>
    <inputParameter>
      <name>eta</name>
      <source>parameters</source>
      <defaultValue>0.3d0</defaultValue>
      <description>The parameter $\eta$ in the \cite{penarrubia_impact_2010} tidal track fitting function.</description>
    </inputParameter>
    <objectBuilder class="darkMatterProfileDMO"  name="darkMatterProfileDMO_"        source="parameters"                                               />
    <objectBuilder class="darkMatterProfileDMO"  name="darkMatterProfileDMOUnheated" source="parameters"   parameterName="darkMatterProfileDMOUnheated"/>
    !!]
    self=outputAnalysisTidalTracksVelocityMaximum(mu,eta,darkMatterProfileDMO_,darkMatterProfileDMOUnheated)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileDMO_"       />
    <objectDestructor name="darkMatterProfileDMOUnheated"/>
    !!]
    return
  end function tidalTracksVelocityMaximumConstructorParameters
  
  function tidalTracksVelocityMaximumConstructorInternal(mu,eta,darkMatterProfileDMO_,darkMatterProfileDMOUnheated) result (self)
    !!{
    Constructor for the \refClass{outputAnalysisTidalTracksVelocityMaximum} output analysis class for internal use.
    !!}
    implicit none
    type            (outputAnalysisTidalTracksVelocityMaximum)                         :: self
    class           (darkMatterProfileDMOClass               ), intent(in   ) , target :: darkMatterProfileDMO_, darkMatterProfileDMOUnheated
    double precision                                          , intent(in   )          :: mu                   , eta
    !![
    <constructorAssign variables="mu, eta, *darkMatterProfileDMO_, *darkMatterProfileDMOUnheated"/>
    !!]

    !$ call OMP_Init_Lock(self%accumulateLock)
    self%countPointsOnTrack=0
    self%logLikelihood_    =0.0d0
    return
  end function tidalTracksVelocityMaximumConstructorInternal

  subroutine tidalTracksVelocityMaximumDestructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisTidalTracksVelocityMaximum} output analysis class.
    !!}
    !$ use :: OMP_Lib, only : OMP_Destroy_Lock
    implicit none
    type(outputAnalysisTidalTracksVelocityMaximum), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileDMO_"       />
    <objectDestructor name="self%darkMatterProfileDMOUnheated"/>
    !!]
    !$ call OMP_Destroy_Lock(self%accumulateLock)
    return
  end subroutine tidalTracksVelocityMaximumDestructor

  subroutine tidalTracksVelocityMaximumAnalyze(self,node,iOutput)
    !!{
    Analyze the maximum velocity tidal track.
    !!}
    use    :: Galacticus_Nodes  , only : nodeComponentBasic   , nodeComponentSatellite
    use    :: Mass_Distributions, only : massDistributionClass
    !$ use :: OMP_Lib           , only : OMP_Set_Lock         , OMP_Unset_Lock
    implicit none
    class           (outputAnalysisTidalTracksVelocityMaximum), intent(inout)             :: self
    type            (treeNode                                ), intent(inout)             :: node
    integer         (c_size_t                                ), intent(in   )             :: iOutput
    class           (nodeComponentBasic                      ), pointer                   :: basic
    class           (nodeComponentSatellite                  ), pointer                   :: satellite
    class           (massDistributionClass                   ), pointer                   :: massDistribution_             , massDistributionUnheated
    double precision                                          , dimension(:), allocatable :: fractionMassBound_            , fractionVelocityMaximum_             , &
         &                                                                                   fractionVelocityMaximumTarget_
    double precision                                                                      :: fractionMassBound             , fractionVelocityMaximum              , &
         &                                                                                   fractionVelocityMaximumTarget , varianceFractionVelocityMaximumTarget

    ! Skip non-satellites.
    if (.not.node%isSatellite()) return
    ! Extract the bound mass and maximum velocity fractions.
    basic                    => node             %basic                           (    )
    satellite                => node             %satellite                       (    )
    massDistribution_        => self             %darkMatterProfileDMO_       %get(node)
    massDistributionUnheated => self             %darkMatterProfileDMOUnheated%get(node)
    fractionMassBound        =  satellite        %boundMass                       (    )/basic                   %mass                        ()
    fractionVelocityMaximum  =  massDistribution_%velocityRotationCurveMaximum    (    )/massDistributionUnheated%velocityRotationCurveMaximum()
    !![
    <objectDestructor name="massDistribution_"       />
    <objectDestructor name="massDistributionUnheated"/>
    !!]
    ! Evaluate the target value. Uses the Penarrubia et al. (2010) fitting function.
    fractionVelocityMaximumTarget        =+ 2.0d0                   **self%mu  &
         &                                *       fractionMassBound **self%eta &
         &                                /(1.0d0+fractionMassBound)**self%mu
    varianceFractionVelocityMaximumTarget=(0.1d0*fractionVelocityMaximumTarget)**2
    !$ call OMP_Set_Lock(self%accumulateLock)
    self%countPointsOnTrack=self%countPointsOnTrack+1
    if (allocated(self%fractionMassBound) .and. size(self%fractionMassBound) < self%countPointsOnTrack) then
       call move_alloc(self%fractionMassBound,            fractionMassBound_            )
       call move_alloc(self%fractionVelocityMaximum      ,fractionVelocityMaximum_      )
       call move_alloc(self%fractionVelocityMaximumTarget,fractionVelocityMaximumTarget_)
       allocate(self%fractionMassBound            (size(fractionMassBound_            )*2))
       allocate(self%fractionVelocityMaximum      (size(fractionVelocityMaximum_      )*2))
       allocate(self%fractionVelocityMaximumTarget(size(fractionVelocityMaximumTarget_)*2))
       self%fractionMassBound            (1:size(fractionMassBound_            ))=fractionMassBound_
       self%fractionVelocityMaximum      (1:size(fractionVelocityMaximum_      ))=fractionVelocityMaximum_
       self%fractionVelocityMaximumTarget(1:size(fractionVelocityMaximumTarget_))=fractionVelocityMaximumTarget_
       deallocate(fractionMassBound_            )
       deallocate(fractionVelocityMaximum_      )
       deallocate(fractionVelocityMaximumTarget_)
    else if (.not.allocated(self%fractionMassBound)) then
       allocate(self%fractionMassBound            (1))
       allocate(self%fractionVelocityMaximum      (1))
       allocate(self%fractionVelocityMaximumTarget(1))
    end if
    self%fractionMassBound            (self%countPointsOnTrack)=fractionMassBound
    self%fractionVelocityMaximum      (self%countPointsOnTrack)=fractionVelocityMaximum
    self%fractionVelocityMaximumTarget(self%countPointsOnTrack)=fractionVelocityMaximumTarget
    self%logLikelihood_                                        =+self%logLikelihood_                   &
         &                                                      -0.5d0                                 &
         &                                                      *(                                     &
         &                                                        +fractionVelocityMaximum             &
         &                                                        -fractionVelocityMaximumTarget       &
         &                                                       )**2                                  &
         &                                                      /varianceFractionVelocityMaximumTarget
    !$ call OMP_Unset_Lock(self%accumulateLock)
    return
  end subroutine tidalTracksVelocityMaximumAnalyze

  subroutine tidalTracksVelocityMaximumReduce(self,reduced)
    !!{
    Reduce over the maximum velocity tidal track output analysis.
    !!}
    use    :: Error  , only : Error_Report
    !$ use :: OMP_Lib, only : OMP_Set_Lock, OMP_Unset_Lock
    implicit none
    class           (outputAnalysisTidalTracksVelocityMaximum), intent(inout)              :: self
    class           (outputAnalysisClass                     ), intent(inout)              :: reduced
    double precision                                          , dimension(:) , allocatable :: fractionMassBound_            , fractionVelocityMaximum_ , &
         &                                                                                    fractionVelocityMaximumTarget_

    select type (reduced)
    class is (outputAnalysisTidalTracksVelocityMaximum)
       !$ call OMP_Set_Lock(reduced%accumulateLock)       
       allocate(fractionMassBound_            (self%countPointsOnTrack+reduced%countPointsOnTrack))
       allocate(fractionVelocityMaximum_      (self%countPointsOnTrack+reduced%countPointsOnTrack))
       allocate(fractionVelocityMaximumTarget_(self%countPointsOnTrack+reduced%countPointsOnTrack))
       if (self%countPointsOnTrack > 0) then
          fractionMassBound_            (                        1:self%countPointsOnTrack             )=self   %fractionMassBound
          fractionVelocityMaximum_      (                        1:self%countPointsOnTrack             )=self   %fractionVelocityMaximum
          fractionVelocityMaximumTarget_(                        1:self%countPointsOnTrack             )=self   %fractionVelocityMaximumTarget
       end if
       if (reduced%countPointsOnTrack > 0) then
          fractionMassBound_            (self%countPointsOnTrack+1:size(fractionMassBound_            ))=reduced%fractionMassBound
          fractionVelocityMaximum_      (self%countPointsOnTrack+1:size(fractionVelocityMaximum_      ))=reduced%fractionVelocityMaximum
          fractionVelocityMaximumTarget_(self%countPointsOnTrack+1:size(fractionVelocityMaximumTarget_))=reduced%fractionVelocityMaximumTarget
          deallocate(reduced%fractionMassBound            )
          deallocate(reduced%fractionVelocityMaximum      )
          deallocate(reduced%fractionVelocityMaximumTarget)
       end if
       call move_alloc(fractionMassBound_            ,reduced%fractionMassBound            )
       call move_alloc(fractionVelocityMaximum_      ,reduced%fractionVelocityMaximum      )
       call move_alloc(fractionVelocityMaximumTarget_,reduced%fractionVelocityMaximumTarget)
       reduced%logLikelihood_    =reduced%logLikelihood_    +self%logLikelihood_
       reduced%countPointsOnTrack=reduced%countPointsOnTrack+self%countPointsOnTrack
       !$ call OMP_Unset_Lock(reduced%accumulateLock)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine tidalTracksVelocityMaximumReduce

  subroutine tidalTracksVelocityMaximumFinalize(self,groupName)
    !!{
    Output results of the maximum velocity tidal track output analysis.
    !!}
#ifdef USEMPI
    use :: MPI_Utilities, only : mpiSelf
#endif
    use :: Output_HDF5  , only : outputFile
    use :: HDF5_Access  , only : hdf5Access
    use :: IO_HDF5      , only : hdf5Object
    implicit none
    class           (outputAnalysisTidalTracksVelocityMaximum), intent(inout)              :: self
    type            (varying_string                          ), intent(in   ), optional    :: groupName
    type            (hdf5Object                              )               , target      :: analysesGroup                 , subGroup
    type            (hdf5Object                              )               , pointer     :: inGroup
    type            (hdf5Object                              )                             :: analysisGroup                 , dataset
#ifdef USEMPI
    integer                                                                                :: offset
    integer                                                   , dimension(:) , allocatable :: countPointsOnTrack
    double precision                                          , dimension(:) , allocatable :: fractionMassBound_            , fractionVelocityMaximum_, &
         &                                                                                    fractionVelocityMaximumTarget_
#endif

#ifdef USEMPI
    ! If running under MPI, accumulate tracks across all processes.
    allocate(countPointsOnTrack(0:mpiSelf%count()-1))
    countPointsOnTrack                =0
    countPointsOnTrack(mpiSelf%rank())=self   %    countPointsOnTrack
    countPointsOnTrack                =mpiSelf%sum(countPointsOnTrack)
    allocate(fractionMassBound_            (sum(countPointsOnTrack)))
    allocate(fractionVelocityMaximum_      (sum(countPointsOnTrack)))
    allocate(fractionVelocityMaximumTarget_(sum(countPointsOnTrack)))
    fractionMassBound_            =0.0d0
    fractionVelocityMaximum_      =0.0d0
    fractionVelocityMaximumTarget_=0.0d0
    if (mpiSelf%rank() == 0) then
       offset=1
    else
       offset=1+sum(countPointsOnTrack(0:mpiSelf%rank()-1))
    end if
    if (self%countPointsOnTrack > 0) then
       fractionMassBound_            (offset:offset+self%countPointsOnTrack-1)=self%fractionMassBound
       fractionVelocityMaximum_      (offset:offset+self%countPointsOnTrack-1)=self%fractionVelocityMaximum
       fractionVelocityMaximumTarget_(offset:offset+self%countPointsOnTrack-1)=self%fractionVelocityMaximumTarget
    end if
    fractionMassBound_            =mpiSelf%sum(fractionMassBound_            )
    fractionVelocityMaximum_      =mpiSelf%sum(fractionVelocityMaximum_      )
    fractionVelocityMaximumTarget_=mpiSelf%sum(fractionVelocityMaximumTarget_)
    if (self%countPointsOnTrack > 0) then
       deallocate(self%fractionMassBound            )
       deallocate(self%fractionVelocityMaximum      )
       deallocate(self%fractionVelocityMaximumTarget)
    end if
    call move_alloc(fractionMassBound_            ,self%fractionMassBound            )
    call move_alloc(fractionVelocityMaximum_      ,self%fractionVelocityMaximum      )
    call move_alloc(fractionVelocityMaximumTarget_,self%fractionVelocityMaximumTarget)
    self%logLikelihood_    =mpiSelf%sum(self%logLikelihood_    )
    self%countPointsOnTrack=        sum(     countPointsOnTrack)
#endif
    !$ call hdf5Access%set()
    analysesGroup =  outputFile   %openGroup('analyses'     )
    inGroup       => analysesGroup
    if (present(groupName)) then
       subGroup   =  analysesGroup%openGroup(char(groupName))
       inGroup    => subGroup
    end if
    analysisGroup=inGroup%openGroup('tidalTrackVelocityMaximum','Analysis of the peak rotation curve velocity vs. bound mass fraction tidal track.')
    ! Write metadata describing this analysis.
    call analysisGroup%writeAttribute('Subhalo peak rotation curve velocity vs. bound mass fraction tidal track.','description'                          )
    call analysisGroup%writeAttribute("function1D"                                                               ,'type'                                 )
    call analysisGroup%writeAttribute('$M_\mathrm{bound}/\mathrm{M}_\odot$'                                      ,'xAxisLabel'                           )
    call analysisGroup%writeAttribute('$V_\mathrm{max}/V_\mathrm{max,0}$'                                        ,'yAxisLabel'                           )
    call analysisGroup%writeAttribute(.true.                                                                     ,'xAxisIsLog'                           )
    call analysisGroup%writeAttribute(.true.                                                                     ,'yAxisIsLog'                           )
    call analysisGroup%writeAttribute('fractionMassBound'                                                        ,'xDataset'                             )
    call analysisGroup%writeAttribute('fractionVelocityMaximum'                                                  ,'yDataset'                             )
    call analysisGroup%writeAttribute('fractionVelocityMaximumTarget'                                            ,'yDatasetTarget'                       )
    call analysisGroup%writeAttribute(self%logLikelihood()                                                       ,'logLikelihood'                        )
    call analysisGroup%writeAttribute('Pe\\~narrubia et al. (2010)'                                              ,'targetLabel'                          )
    call analysisGroup%writeDataset  (self%fractionMassBound            ,'fractionMassBound'            ,'Fraction of bound mass remaining'           ,datasetReturned=dataset)
    call dataset      %writeAttribute(' '                               ,'units'                                                                                              )
    call dataset      %writeAttribute(1.0d0                             ,'unitsInSI'                                                                                          )
    call analysisGroup%writeDataset  (self%fractionVelocityMaximum      ,'fractionVelocityMaximum'      ,'Fraction of peak rotation velocity'         ,datasetReturned=dataset)
    call dataset      %writeAttribute(' '                               ,'units'                                                                                              )
    call dataset      %writeAttribute(1.0d0                             ,'unitsInSI'                                                                                          )
    call analysisGroup%writeDataset  (self%fractionVelocityMaximumTarget,'fractionVelocityMaximumTarget','Fraction of peak rotation velocity [target]',datasetReturned=dataset)
    call dataset      %writeAttribute(' '                               ,'units'                                                                                              )
    call dataset      %writeAttribute(1.0d0                             ,'unitsInSI'                                                                                          )
    !$ call hdf5Access%unset()
    return
  end subroutine tidalTracksVelocityMaximumFinalize
  
  double precision function tidalTracksVelocityMaximumLogLikelihood(self)
    !!{
    Return the log-likelihood of a maximum velocity tidal track output analysis.
    !!}
    implicit none
    class(outputAnalysisTidalTracksVelocityMaximum), intent(inout) :: self

    tidalTracksVelocityMaximumLogLikelihood=self%logLikelihood_
    return
  end function tidalTracksVelocityMaximumLogLikelihood

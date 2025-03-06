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
  Implements an output analysis class that evaluates model likelihood given the galaxy mass-size relations of \cite{shen_size_2003}.
  !!}

  !$ use :: OMP_Lib            , only : omp_lock_kind
  use    :: Cosmology_Functions, only : cosmologyFunctionsClass
  
  !![
  <outputAnalysis name="outputAnalysisMassSizeRelationShen2003">
    <description>An output analysis class that evaluates model likelihood given the galaxy mass-size relations of \cite{shen_size_2003}.</description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisClass) :: outputAnalysisMassSizeRelationShen2003
     !!{
     An output analysis class that evaluates model likelihood given the galaxy mass-size relations of \cite{shen_size_2003}.
     !!}
     private
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null()
     !$ integer      (omp_lock_kind          )          :: accumulateLock
     double precision                                   :: timeMinimum                  , timeMaximum, &
          &                                                logLikelihood_
     logical                                            :: finalized
   contains
     !![
     <methods>
       <method method="finalizeAnalysis" description="Finalize the analysis of this function."/>
     </methods>
     !!]
     final     ::                     massSizeRelationShen2003Destructor
     procedure :: analyze          => massSizeRelationShen2003Analyze
     procedure :: finalize         => massSizeRelationShen2003Finalize
     procedure :: finalizeAnalysis => massSizeRelationShen2003FinalizeAnalysis
     procedure :: reduce           => massSizeRelationShen2003Reduce
     procedure :: logLikelihood    => massSizeRelationShen2003LogLikelihood
  end type outputAnalysisMassSizeRelationShen2003

  interface outputAnalysisMassSizeRelationShen2003
     !!{
     Constructors for the ``massSizeRelationShen2003'' output analysis class.
     !!}
     module procedure massSizeRelationShen2003ConstructorParameters
     module procedure massSizeRelationShen2003ConstructorInternal
  end interface outputAnalysisMassSizeRelationShen2003

contains

  function massSizeRelationShen2003ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``massSizeRelationShen2003'' output analysis class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (outputAnalysisMassSizeRelationShen2003)                :: self
    type (inputParameters                       ), intent(inout) :: parameters
    class(cosmologyFunctionsClass               ), pointer       :: cosmologyFunctions_ => null()
    
    !![
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    self=outputAnalysisMassSizeRelationShen2003(cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function massSizeRelationShen2003ConstructorParameters
  
  function massSizeRelationShen2003ConstructorInternal(cosmologyFunctions_) result (self)
    !!{
    Constructor for the ``massSizeRelationShen2003'' output analysis class for internal use.
    !!}
    implicit none
    type (outputAnalysisMassSizeRelationShen2003)                        :: self
    class(cosmologyFunctionsClass               ), intent(in   ), target :: cosmologyFunctions_
    !![
    <constructorAssign variables="*cosmologyFunctions_"/>
    !!]
    
    !$ call OMP_Init_Lock(self%accumulateLock)
    self%logLikelihood_=0.0d0
    self%finalized     =.false.
    ! Find the span of time for the SDSS sample - this is a crude selection motivated by the span of redshifts in Figure 1 of Shen
    ! et al. (2003).
    self%timeMinimum=self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(0.3d0))
    self%timeMaximum=self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(0.0d0))
    return
  end function massSizeRelationShen2003ConstructorInternal

  subroutine massSizeRelationShen2003Destructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily massSizeRelationShen2003} output analysis class.
    !!}
    !$ use :: OMP_Lib, only : OMP_Destroy_Lock
    implicit none
    type(outputAnalysisMassSizeRelationShen2003), intent(inout) :: self

    !$ call OMP_Destroy_Lock(self%accumulateLock)
    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine massSizeRelationShen2003Destructor

  subroutine massSizeRelationShen2003Analyze(self,node,iOutput)
    !!{
    Analyze the maximum velocity tidal track.
    !!}
    use    :: Galacticus_Nodes            , only : nodeComponentBasic
    use    :: Galactic_Structure_Options  , only : componentTypeDisk    , componentTypeSpheroid, massTypeStellar
    use    :: Mass_Distributions          , only : massDistributionClass
    use    :: Models_Likelihoods_Constants, only : logImpossible
    use    :: Numerical_Constants_Math    , only : Pi
    use    :: Numerical_Constants_Prefixes, only : mega                 , kilo
    !$ use :: OMP_Lib                     , only : OMP_Set_Lock         , OMP_Unset_Lock
    implicit none
    class           (outputAnalysisMassSizeRelationShen2003), intent(inout) :: self
    type            (treeNode                              ), intent(inout) :: node
    integer         (c_size_t                              ), intent(in   ) :: iOutput
    class           (massDistributionClass                 ), pointer       :: massDistributionDisk        , massDistributionSpheroid         , &
         &                                                                     massDistribution_
    class           (nodeComponentBasic                    ), pointer       :: basic
    double precision                                        , parameter     :: a                    =0.56d0, b                       =2.88d-06, &
         &                                                                     alpha                =0.14d0, beta                    =0.39d+00, &
         &                                                                     gamma                =0.10d0, massZeroPoint           =3.98d+10, &
         &                                                                     sigma1               =0.47d0, sigma2                  =0.34d+00
    double precision                                                        :: massDisk                    , massSpheroid                     , &
         &                                                                     radius                      , mass                             , &
         &                                                                     massLogarithmic             , radiusLogarithmic                , &
         &                                                                     radiusLogarithmicMean       , scatterLogarithmic

    ! Skip satellites.
    if (node%isSatellite()) return
    ! Skip galaxies outside of the relevant redshift interval.
    basic => node%basic()
    if     (                                 &
         &   basic%time() < self%timeMinimum &
         &  .or.                             &
         &   basic%time() > self%timeMaximum &
         & ) return
    ! Extract the stellar masses.
    massDistribution_        => node                    %massDistribution   (                                     massType=massTypeStellar)
    massDistributionDisk     => node                    %massDistribution   (componentType =componentTypeDisk    ,massType=massTypeStellar)
    massDistributionSpheroid => node                    %massDistribution   (componentType =componentTypeSpheroid,massType=massTypeStellar)
    massDisk                 =  massDistributionDisk    %massTotal          (                                                             )
    massSpheroid             =  massDistributionSpheroid%massTotal          (                                                             )
    radius                   =  massDistribution_       %radiusEnclosingMass(massFractional=0.5d0                                         )
    mass                     =  +massDisk     &
         &                      +massSpheroid
    !![
    <objectDestructor name="massDistribution_"       />
    <objectDestructor name="massDistributionDisk"    />
    <objectDestructor name="massDistributionSpheroid"/>
    !!]
    !$ call OMP_Set_Lock(self%accumulateLock)
    if     (                                     &
         &   self%logLikelihood_ > logImpossible &
         &  .and.                                &
         &        mass           > 0.0d0         &
         &  .and.                                &
         &        radius         > 0.0d0         &
         &) then
       ! Get logarithmic mass and size.
       massLogarithmic  =log10(mass  )
       radiusLogarithmic=log  (radius)
       ! Evaluate the parameters of the Shen et al. (2003) distribution function. (Convert from kpc to Mpc.)
       if (massDisk > 0.5d0*mass) then
          ! Late-type galaxy (Shen et al. 2003; eq. 18).
          radiusLogarithmicMean=log(gamma)+alpha*log(mass)+(beta-alpha)*log(1.0d0+mass/massZeroPoint)+log(kilo/mega)
       else
          ! Early-type galaxy (Shen et al. 2003; eq. 17).
          radiusLogarithmicMean=log(    b)+    a*log(mass)                                           +log(kilo/mega)
       end if
       scatterLogarithmic   =sigma2+(sigma1-sigma2)/(1.0d0+(mass/massZeroPoint)**2) ! Shen et al. (2003; eq. 19).
       ! Evaluate the likelihood.
       self%logLikelihood_=+self%logLikelihood_        &
            &              -0.5d0                      &
            &              *   (                       &
            &                   +radiusLogarithmic     &
            &                   -radiusLogarithmicMean &
            &                  )**2                    &
            &              /     scatterLogarithmic**2 &
            &              -0.5d0                      &
            &              *log(                       &
            &                   +2.0d0                 &
            &                   *Pi                    &
            &                   *scatterLogarithmic**2 &
            &                  )
    else
       self%logLikelihood_=logImpossible
    end if
    !$ call OMP_Unset_Lock(self%accumulateLock)
    return
  end subroutine massSizeRelationShen2003Analyze

  subroutine massSizeRelationShen2003Reduce(self,reduced)
    !!{
    Reduce over the mass-size output analysis.
    !!}
    use    :: Error  , only : Error_Report
    !$ use :: OMP_Lib, only : OMP_Set_Lock, OMP_Unset_Lock
    implicit none
    class(outputAnalysisMassSizeRelationShen2003), intent(inout) :: self
    class(outputAnalysisClass                   ), intent(inout) :: reduced

    select type (reduced)
    class is (outputAnalysisMassSizeRelationShen2003)
       !$ call OMP_Set_Lock(reduced%accumulateLock)       
       reduced%logLikelihood_=reduced%logLikelihood_+self%logLikelihood_
       !$ call OMP_Unset_Lock(reduced%accumulateLock)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine massSizeRelationShen2003Reduce

  subroutine massSizeRelationShen2003FinalizeAnalysis(self)
    !!{
    Compute final likelihood.
    !!}
#ifdef USEMPI
    use :: MPI_Utilities, only : mpiSelf
#endif
    implicit none
    class(outputAnalysisMassSizeRelationShen2003), intent(inout) :: self

    ! If already finalized, no need to do anything.
    if (self%finalized) return
    self%finalized=.true.
#ifdef USEMPI
    ! If running under MPI, accumulate likelihood across all processes.
    self%logLikelihood_=mpiSelf%sum(self%logLikelihood_)
#endif
    return
  end subroutine massSizeRelationShen2003FinalizeAnalysis

  subroutine massSizeRelationShen2003Finalize(self,groupName)
    !!{
    Output results of the mass-size relation output analysis.
    !!}
    use :: Output_HDF5  , only : outputFile
    use :: HDF5_Access  , only : hdf5Access
    use :: IO_HDF5      , only : hdf5Object
    implicit none
    class(outputAnalysisMassSizeRelationShen2003), intent(inout)           :: self
    type (varying_string                        ), intent(in   ), optional :: groupName
    type (hdf5Object                            )               , target   :: analysesGroup, subGroup
    type (hdf5Object                            )               , pointer  :: inGroup
    type (hdf5Object                            )                          :: analysisGroup

    ! Finalize analysis.
    call self%finalizeAnalysis()
    !$ call hdf5Access%set()
    analysesGroup =  outputFile   %openGroup('analyses'     )
    inGroup       => analysesGroup
    if (present(groupName)) then
       subGroup   =  analysesGroup%openGroup(char(groupName))
       inGroup    => subGroup
    end if
    analysisGroup=inGroup%openGroup('massSizeRelationShen2003','Analysis of the mass-size relation of Shen et al. (2003).')
    ! Write metadata describing this analysis.
    call analysisGroup%writeAttribute('The mass-size relation of Shen et al. (2003).','description'  )
    call analysisGroup%writeAttribute(self%logLikelihood()                           ,'logLikelihood')
    call analysisGroup%close         (                                                               )
    if (present(groupName)) &
         & call subGroup %close      (                                                               )
    call analysesGroup%close         (                                                               )
    !$ call hdf5Access%unset()
    return
  end subroutine massSizeRelationShen2003Finalize
  
  double precision function massSizeRelationShen2003LogLikelihood(self)
    !!{
    Return the log-likelihood of the mass-size output analysis.
    !!}
    implicit none
    class(outputAnalysisMassSizeRelationShen2003), intent(inout) :: self

    call self%finalizeAnalysis()
    massSizeRelationShen2003LogLikelihood=self%logLikelihood_
    return
  end function massSizeRelationShen2003LogLikelihood

!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

!% Defines methods for a mass function likelihood function.

function likelihoodHaloSpinDistributionConstructor(fileName,distributionName,redshift,massHaloMinimum,massParticle,particleCountMinimum,energyEstimateParticleCountMaximum) result (self)
  !% Constructor for the halo mass function likelihood class.
  use IO_HDF5
  use Memory_Management
  use Cosmology_Functions
  use Galacticus_Error
  implicit none
  type            (likelihoodHaloSpinDistribution)                :: self
  character       (len=*                         ), intent(in   ) :: fileName               , distributionName
  double precision                                , intent(in   ) :: redshift               , massHaloMinimum                   , &
       &                                                             massParticle           , energyEstimateParticleCountMaximum
  integer                                         , intent(in   ) :: particleCountMinimum
  class           (cosmologyFunctionsClass       ), pointer       :: cosmologyFunctions_
  type            (hdf5Object                    )                :: spinDistributionFile
  double precision                                                :: spinIntervalLogarithmic
  integer                                                         :: i
  
  ! Convert redshift to time.
  cosmologyFunctions_ => cosmologyFunctions()
  self%time           =  cosmologyFunctions_ %cosmicTime                 (          &
       &                  cosmologyFunctions_%expansionFactorFromRedshift (         &
       &                                                                   redshift &
       &                                                                  )         &
       &                                                                 )
  ! Store particle mass and minimum halo mass.
  self%massParticle                      =massParticle
  self%massHaloMinimum                   =massHaloMinimum
  self%particleCountMinimum              =particleCountMinimum
  self%energyEstimateParticleCountMaximum=energyEstimateParticleCountMaximum
  ! Determine the type of spin distribution.
  select case (trim(distributionName))
  case ("logNormal")
     self%distributionType=distributionTypeLogNormal
  case ("bett2007")
     self%distributionType=distributionTypeBett2007
  case default
     call Galacticus_Error_Report('likelihoodHaloSpinDistributionConstructor','unknown distribution type')
  end select
  ! Read the target spin distribution from file.
  !$omp critical(HDF5_Access)
  call spinDistributionFile%openFile   (trim(fileName),readOnly=.true.)
  call spinDistributionFile%readDataset("spinParameter"    ,self%spin             )
  call spinDistributionFile%readDataset("distribution"     ,self%distribution     )
  call spinDistributionFile%readDataset("distributionError",self%distributionError)
  call spinDistributionFile%close()
  !$omp end critical(HDF5_Access)
  ! Compute spin ranges for bins.
  spinIntervalLogarithmic=+log(                            &
       &                       +self%spin(size(self%spin)) &
       &                       /self%spin(              1) &
       &                      )                            &
       &                  /dble(                           &
       &                        +size(self%spin)           &
       &                        -1                         &
       &                       )
  call allocateArray(self%spinMinimum,shape(self%spin))
  call allocateArray(self%spinMaximum,shape(self%spin))
  do i=1,size(self%spin)
     self%spinMinimum(i)=self%spin(i)*exp(-0.5d0*spinIntervalLogarithmic)
     self%spinMaximum(i)=self%spin(i)*exp(+0.5d0*spinIntervalLogarithmic)
  end do
  return
end function likelihoodHaloSpinDistributionConstructor

double precision function likelihoodHaloSpinDistributionEvaluate(self,simulationState,parameterMappings,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed,timeEvaluate,logLikelihoodVariance)
  !% Return the log-likelihood for the halo spin distribution likelihood function.
  use, intrinsic :: ISO_C_Binding
  use               Numerical_Integration
  use               MPI_Utilities
  use               Constraints_State
  use               Constraints_Constants
  use               Constraints_Convergence
  use               Galacticus_Error
  use               Halo_Mass_Functions
  use               Statistics_NBody_Halo_Mass_Errors
  use               Halo_Spin_Distributions
  use               Galacticus_Nodes
  use               Dark_Matter_Profiles
  use               Dark_Matter_Halo_Scales
  implicit none
  class           (likelihoodHaloSpinDistribution ), intent(inout)               :: self
  class           (state                          ), intent(inout)               :: simulationState
  type            (mappingList                    ), intent(in   ), dimension(:) :: parameterMappings
  class           (convergence                    ), intent(inout)               :: simulationConvergence
  double precision                                 , intent(in   )               :: temperature          , logLikelihoodCurrent, &
       &                                                                            logPriorCurrent      , logPriorProposed
  real                                             , intent(inout)               :: timeEvaluate
  double precision                                 , intent(  out), optional     :: logLikelihoodVariance
  double precision                                 , allocatable  , dimension(:) :: stateVector          , distribution
  class           (haloMassFunctionClass          ), pointer                     :: haloMassFunction_
  class           (nbodyHaloMassErrorClass        ), pointer                     :: nbodyHaloMassError_
  class           (darkMatterProfileClass         ), pointer                     :: darkMatterProfile_
  class           (darkMatterHaloScaleClass       ), pointer                     :: darkMatterHaloScale_
  type            (treeNode                       ), pointer                     :: node
  class           (nodeComponentBasic             ), pointer                     :: nodeBasic
  class           (nodeComponentSpin              ), pointer                     :: nodeSpin
  class           (haloSpinDistributionLogNormal  ), pointer                     :: distributionLogNormal
  class           (haloSpinDistributionBett2007   ), pointer                     :: distributionBett2007
  type            (haloSpinDistributionNbodyErrors)                              :: distributionNbody
  type            (fgsl_function                  )                              :: integrandFunction
  type            (fgsl_integration_workspace     )                              :: integrationWorkspace
  integer                                                                        :: i
  !GCC$ attributes unused :: simulationConvergence, temperature, timeEvaluate, logLikelihoodCurrent, logPriorCurrent
  
  ! There is no variance in our likelihood estimate.
  if (present(logLikelihoodVariance)) logLikelihoodVariance=0.0d0
  ! Do not evaluate if the proposed prior is impossible.
  if (logPriorProposed <= logImpossible) then
     likelihoodHaloSpinDistributionEvaluate=0.0d0
     return
  end if
  ! Get required objects.
  haloMassFunction_    => haloMassFunction   ()
  nbodyHaloMassError_  => nbodyHaloMassError ()
  darkMatterHaloScale_ => darkMatterHaloScale()
  darkMatterProfile_   => darkMatterProfile  ()
  ! Extract parameters of the spin distribution.
  stateVector=simulationState%get()
  do i=1,size(stateVector)
     stateVector(i)=parameterMappings(i)%thisMapping%unmap(stateVector(i))
  end do
  ! Build the halo spin distribution object.
  select case (self%distributionType)
  case (distributionTypeLogNormal)
     if (size(stateVector) /= 2)                                                                            &
          & call Galacticus_Error_Report(                                                                   & 
          &                              'likelihoodHaloSpinDistributionEvaluate'                         , &
          &                              '2 parameters are required for the "lognormal" spin distribution'  &
          &                             )
     allocate(haloSpinDistributionLogNormal :: distributionLogNormal)
     select type (distributionLogNormal)
     type is (haloSpinDistributionLogNormal)
        distributionLogNormal=haloSpinDistributionLogNormal  (                                             &
             &                                                stateVector(1)                             , &
             &                                                stateVector(2)                               &
             &                                               )
        distributionNbody    =haloSpinDistributionNbodyErrors(                                             &
             &                                                distributionLogNormal                      , &
             &                                                self%massParticle                          , &
             &                                                self%particleCountMinimum                  , &
             &                                                self%energyEstimateParticleCountMaximum    , &
             &                                                self%time                                  , &
             &                                                nbodyHaloMassError_                        , &
             &                                                haloMassFunction_                          , &
             &                                                darkMatterHaloScale_                       , &
             &                                                darkMatterProfile_                           &
             &                                               )
     end select
  case (distributionTypeBett2007)
     if (size(stateVector) /= 2)                                                                           &
          & call Galacticus_Error_Report(                                                                  & 
          &                              'likelihoodHaloSpinDistributionEvaluate'                        , &
          &                              '2 parameters are required for the "bett2007" spin distribution'  &
          &                             )
     allocate(haloSpinDistributionBett2007 :: distributionBett2007)
     select type (distributionBett2007)
     type is (haloSpinDistributionBett2007)
        distributionBett2007=haloSpinDistributionBett2007    (                                             &
             &                                                stateVector(1)                             , &
             &                                                stateVector(2)                               &
             &                                               )
        distributionNbody    =haloSpinDistributionNbodyErrors(                                             &
             &                                                distributionBett2007                       , &
             &                                                self%massParticle                          , &
             &                                                self%particleCountMinimum                  , &
             &                                                self%energyEstimateParticleCountMaximum    , &
             &                                                self%time                                  , &
             &                                                nbodyHaloMassError_                        , &
             &                                                haloMassFunction_                          , &
             &                                                darkMatterHaloScale_                       , &
             &                                                darkMatterProfile_                           &
             &                                               )
     end select
  end select
  ! Nullify objects to avoid them being destroyed.
  nullify(haloMassFunction_    )
  nullify(nbodyHaloMassError_  )
  nullify(darkMatterHaloScale_ )
  nullify(darkMatterProfile_   )
  nullify(distributionLogNormal)
  ! Create a work node and set the appropriate cosmological time.
  node      => treeNode      (                 )
  nodeBasic => node    %basic(autoCreate=.true.)
  nodeSpin  => node    %spin (autoCreate=.true.)
  call nodeBasic%timeSet(self%time)
  ! Compute the spin distribution. The spin distribution is averaged over the width of each bin.
  allocate(distribution(size(self%spin)))
  do i=1,size(self%spin)
     distribution(i)=+Integrate(                                  &
          &                     self%spinMinimum(i)             , &
          &                     self%spinMaximum(i)             , &
          &                     spinDistributionIntegrate       , &
          &                     integrandFunction               , &
          &                     integrationWorkspace            , &
          &                     toleranceAbsolute        =0.0d+0, &
          &                     toleranceRelative        =1.0d-6  &
          &                    )                                  &
          &          /log10(                                      &
          &                 +self%spinMaximum(i)                  &
          &                 /self%spinMinimum(i)                  &
          &                )
     call Integrate_Done(integrandFunction,integrationWorkspace)
  end do
  ! Evaluate the log-likelihood.
  likelihoodHaloSpinDistributionEvaluate=-0.5d0*sum(((distribution-self%distribution)/self%distributionError)**2,mask=self%distributionError > 0.0d0)
  ! Clean up.
  call node%destroy()
  deallocate(node        )
  deallocate(stateVector )
  deallocate(distribution)  
  return

contains

  double precision function spinDistributionIntegrate(spinPrime)
    !% Integrand function used to find cumulative spin distribution over a bin.
    implicit none
    double precision, intent(in   ) :: spinPrime

    call nodeSpin%spinSet(spinPrime)
    spinDistributionIntegrate=distributionNbody%distributionAveraged(node,self%massHaloMinimum)
    return
  end function spinDistributionIntegrate

end function likelihoodHaloSpinDistributionEvaluate

subroutine likelihoodHaloSpinDistributionFunctionChanged(self)
  !% Respond to possible changes in the likelihood function.
  implicit none
  class(likelihoodHaloSpinDistribution), intent(inout) :: self
  !GCC$ attributes unused :: self

  return
end subroutine likelihoodHaloSpinDistributionFunctionChanged

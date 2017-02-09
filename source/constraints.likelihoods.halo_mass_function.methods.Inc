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

function likelihoodHaloMassFunctionConstructor(fileName,redshift,massRangeMinimum,binCountMinimum,massFunctionType,massFunctionErrorModel,massParticle)
  !% Constructor for the halo mass function likelihood class.
  use IO_HDF5
  use Memory_Management
  use Cosmology_Functions
  use Galacticus_Display
  use Galacticus_Error
  implicit none
  type            (likelihoodHaloMassFunction)                                :: likelihoodHaloMassFunctionConstructor
  character       (len=*                     ), intent(in   )                 :: fileName                             , massFunctionType , &
       &                                                                         massFunctionErrorModel
  double precision                            , intent(in   )                 :: redshift                             , massRangeMinimum , &
       &                                                                         massParticle
  integer                                     , intent(in   )                 :: binCountMinimum
  class           (cosmologyFunctionsClass   ), pointer                       :: cosmologyFunctions_
  double precision                            , allocatable  , dimension(:  ) :: eigenValueArray                      , massOriginal     , &
       &                                                                         massFunctionOriginal
  double precision                            , allocatable  , dimension(:,:) :: massFunctionCovarianceOriginal
  character       (len=12                    )                                :: redshiftLabel                        , typeLabel
  type            (hdf5Object                )                                :: massFunctionFile                     , massFunctionGroup, &
       &                                                                         analysisGroup
  integer                                                                     :: i                                    , j                , &
       &                                                                         ii                                   , jj               , &
       &                                                                         massCountReduced
  double precision                                                            :: massIntervalLogarithmic
  type            (matrix                    )                                :: eigenVectors
  type            (vector                    )                                :: eigenValues

  ! Convert redshift to time.
  cosmologyFunctions_ => cosmologyFunctions()
  likelihoodHaloMassFunctionConstructor%time=cosmologyFunctions_ %cosmicTime                 (          &
       &                                      cosmologyFunctions_%expansionFactorFromRedshift (         &
       &                                                                                       redshift &
       &                                                                                      )         &
       &                                                                                     )
  ! Read the halo mass function file.
  write (redshiftLabel,'(f6.3)') redshift
  select case (trim(massFunctionType))
  case ("regular" )
     typeLabel="halo"
  case ("isolated")
     typeLabel="isolatedHalo"
  case default
     call Galacticus_Error_Report('likelihoodHaloMassFunctionConstructor','unknown mass function type')
  end select
  !$omp critical(HDF5_Access)
  call massFunctionFile %openFile(trim(fileName),readOnly=.true.)
  analysisGroup    =massFunctionFile%openGroup('analysis'                                               )
  massFunctionGroup=analysisGroup   %openGroup(trim(typeLabel)//'MassFunctionZ'//trim(adjustl(redshiftLabel)))
  call massFunctionGroup%readDataset("mass"                  ,massOriginal                  )
  call massFunctionGroup%readDataset("massFunction"          ,massFunctionOriginal          )
  call massFunctionGroup%readDataset("massFunctionCovariance",massFunctionCovarianceOriginal)
  call massFunctionGroup%close()
  call analysisGroup    %close()
  call massFunctionFile %close()
  !$omp end critical(HDF5_Access)
  ! Find a reduced mass function excluding any empty bins.
  massCountReduced=0
  do i=1,size(massOriginal)
     if (     massFunctionOriginal          (i  )  <= 0.0d0                                              ) cycle
     if (     massOriginal                  (i  )  <= massRangeMinimum                                   ) cycle
     if (sqrt(massFunctionCovarianceOriginal(i,i)) >= massFunctionOriginal(i)/sqrt(dble(binCountMinimum))) cycle
     massCountReduced=massCountReduced+1
  end do
  if (massCountReduced == 0) call Galacticus_Error_Report('likelihoodHaloMassFunctionConstructor','no usable bins in mass function')  
  call allocateArray(likelihoodHaloMassFunctionConstructor%mass            ,[massCountReduced                 ])
  call allocateArray(likelihoodHaloMassFunctionConstructor%massFunction    ,[massCountReduced                 ])
  call allocateArray(likelihoodHaloMassFunctionConstructor%covarianceMatrix,[massCountReduced,massCountReduced])
  ii=0
  do i=1,size(massOriginal)
     if (     massFunctionOriginal          (i  )  <= 0.0d0                                              ) cycle
     if (     massOriginal                  (i  )  <= massRangeMinimum                                   ) cycle
     if (sqrt(massFunctionCovarianceOriginal(i,i)) >= massFunctionOriginal(i)/sqrt(dble(binCountMinimum))) cycle
     ii=ii+1
     likelihoodHaloMassFunctionConstructor%mass        (ii)=massOriginal        (i)
     likelihoodHaloMassFunctionConstructor%massFunction(ii)=massFunctionOriginal(i)
     jj=0
     do j=1,size(massOriginal)
        if (     massFunctionOriginal          (j  )  <= 0.0d0                                              ) cycle
        if (     massOriginal                  (j  )  <= massRangeMinimum                                   ) cycle
        if (sqrt(massFunctionCovarianceOriginal(j,j)) >= massFunctionOriginal(j)/sqrt(dble(binCountMinimum))) cycle
        jj=jj+1
        likelihoodHaloMassFunctionConstructor%covarianceMatrix(ii,jj)=massFunctionCovarianceOriginal(i,j)        
     end do
  end do
  ! Compute mass ranges for bins.
  massIntervalLogarithmic=+log(                                  &
       &                       +massOriginal(size(massOriginal)) &
       &                       /massOriginal(                 1) &
       &                      )                                  &
       &                  /dble(                                 &
       &                        +size(massOriginal)              &
       &                        -1                               &
       &                       )
  call allocateArray(likelihoodHaloMassFunctionConstructor%massMinimum,shape(likelihoodHaloMassFunctionConstructor%mass))
  call allocateArray(likelihoodHaloMassFunctionConstructor%massMaximum,shape(likelihoodHaloMassFunctionConstructor%mass))
  do i=1,size(likelihoodHaloMassFunctionConstructor%mass)
     likelihoodHaloMassFunctionConstructor%massMinimum(i)=likelihoodHaloMassFunctionConstructor%mass(i)*exp(-0.5d0*massIntervalLogarithmic)
     likelihoodHaloMassFunctionConstructor%massMaximum(i)=likelihoodHaloMassFunctionConstructor%mass(i)*exp(+0.5d0*massIntervalLogarithmic)
  end do
  ! Find the inverse covariance matrix.
  likelihoodHaloMassFunctionConstructor%covariance       =likelihoodHaloMassFunctionConstructor%covarianceMatrix
  likelihoodHaloMassFunctionConstructor%inverseCovariance=likelihoodHaloMassFunctionConstructor%covariance      %invert()
  ! Symmetrize the inverse covariance matrix.
  call likelihoodHaloMassFunctionConstructor%inverseCovariance%symmetrize()
  ! Get eigenvalues and vectors of the inverse covariance matrix.
  allocate(eigenValueArray(size(likelihoodHaloMassFunctionConstructor%mass)))
  call likelihoodHaloMassFunctionConstructor%inverseCovariance%eigenSystem(eigenVectors,eigenValues)
  eigenValueArray=eigenValues
  if (any(eigenValueArray < 0.0d0)) call Galacticus_Display_Message('WARNING: inverse covariance matrix is not semi-positive definite')
  deallocate(eigenValueArray               )
  deallocate(massOriginal                  )
  deallocate(massFunctionOriginal          )
  deallocate(massFunctionCovarianceOriginal)
  ! Determine the error model to use.
  select case (trim(massFunctionErrorModel))
  case ("none"                )
     likelihoodHaloMassFunctionConstructor%errorModel  =haloMassFunctionErrorModelNone
  case ("powerLaw"            )
     likelihoodHaloMassFunctionConstructor%errorModel  =haloMassFunctionErrorModelPowerLaw
  case ("sphericalOverdensity")
     likelihoodHaloMassFunctionConstructor%errorModel  =haloMassFunctionErrorModelSphericalOverdensity
     likelihoodHaloMassFunctionConstructor%massParticle=massParticle
  end select
  return
end function likelihoodHaloMassFunctionConstructor

double precision function likelihoodHaloMassFunctionEvaluate(self,simulationState,parameterMappings,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed,timeEvaluate,logLikelihoodVariance)
  !% Return the log-likelihood for the halo mass function likelihood function.
  use MPI_Utilities
  use Constraints_State
  use Constraints_Constants
  use Constraints_Convergence
  use Galacticus_Error
  use Halo_Mass_Functions
  use Dark_Matter_Halo_Scales
  use Dark_Matter_Profiles
  use Cosmology_Parameters
  use Critical_Overdensities
  use Cosmological_Mass_Variance
  use Statistics_NBody_Halo_Mass_Errors
  implicit none
  class           (likelihoodHaloMassFunction    ), intent(inout)               :: self
  class           (state                         ), intent(in   )               :: simulationState
  type            (mappingList                   ), intent(in   ), dimension(:) :: parameterMappings
  class           (convergence                   ), intent(inout)               :: simulationConvergence
  double precision                                , intent(in   )               :: temperature                     , logLikelihoodCurrent, &
       &                                                                           logPriorCurrent                 , logPriorProposed
  real                                            , intent(inout)               :: timeEvaluate
  double precision                                , intent(  out), optional     :: logLikelihoodVariance
  double precision                                , allocatable  , dimension(:) :: stateVector                     , massFunction
  double precision                                , parameter                   :: errorFractionalMaximum    =0.9d0
  class           (cosmologyParametersClass      ), pointer                     :: cosmologyParameters_
  class           (cosmologicalMassVarianceClass ), pointer                     :: cosmologicalMassVariance_
  class           (criticalOverdensityClass      ), pointer                     :: criticalOverdensity_
  class           (darkMatterHaloScaleClass      ), pointer                     :: darkMatterHaloScale_
  class           (darkMatterProfileClass        ), pointer                     :: darkMatterProfile_
  class           (haloMassFunctionClass         ), pointer                     :: haloMassFunction_
  class           (nbodyHaloMassErrorClass       ), pointer                     :: nbodyHaloMassError_
  type            (haloMassFunctionShethTormen   ), target                      :: haloMassFunctionRaw_
  type            (haloMassFunctionErrorConvolved), target                      :: haloMassFunctionConvolved_
  type            (vector                        )                              :: difference
  integer                                                                       :: i
  !GCC$ attributes unused :: simulationConvergence, temperature, timeEvaluate, logLikelihoodCurrent, logPriorCurrent
  
  ! There is no variance in our likelihood estimate.
  if (present(logLikelihoodVariance)) logLikelihoodVariance=0.0d0
  ! Do not evaluate if the proposed prior is impossible.
  if (logPriorProposed <= logImpossible) then
     likelihoodHaloMassFunctionEvaluate=0.0d0
     return
  end if
  ! Get required objects.
  cosmologyParameters_      => cosmologyParameters     ()
  cosmologicalMassVariance_ => cosmologicalMassVariance()
  criticalOverdensity_      => criticalOverdensity     ()
  ! Build the halo mass function object.
  stateVector=simulationState%get()
  do i=1,size(stateVector)
     stateVector(i)=parameterMappings(i)%thisMapping%unmap(stateVector(i))
  end do
  ! Construct mass function.
  select case (self%errorModel)
  case (haloMassFunctionErrorModelNone                )
     if (size(stateVector) /= 3 )                                                                  &
          & call Galacticus_Error_Report(                                                          & 
          &                              'likelihoodHaloMassFunctionEvaluate'                    , &
          &                              '3 parameters are required for this likelihood function'  &
          &                             )
     ! Use an unconvolved mass function.
     haloMassFunctionRaw_=haloMassFunctionShethTormen(                              &
          &                                           cosmologyParameters_        , &
          &                                           cosmologicalMassVariance_   , &
          &                                           criticalOverdensity_        , &
          &                                           stateVector              (1), &
          &                                           stateVector              (2), &
          &                                           stateVector              (3)  &
          &                                          )
     haloMassFunction_ => haloMassFunctionRaw_
  case (haloMassFunctionErrorModelPowerLaw            )
     if (size(stateVector) /= 5 )                                                                  &
          & call Galacticus_Error_Report(                                                          & 
          &                              'likelihoodHaloMassFunctionEvaluate'                    , &
          &                              '3 parameters are required for this likelihood function'  &
          &                             )
     ! Use a mass function convolved with an error distribution.
     allocate(haloMassFunctionShethTormen::haloMassFunction_)
     select type (haloMassFunction_)
     type is (haloMassFunctionShethTormen)
        haloMassFunction_         =haloMassFunctionShethTormen   (                              &
             &                                                    cosmologyParameters_        , &
             &                                                    cosmologicalMassVariance_   , &
             &                                                    criticalOverdensity_        , &
             &                                                    stateVector              (1), &
             &                                                    stateVector              (2), &
             &                                                    stateVector              (3)  &
             &                                                   )
     end select
     allocate(nbodyHaloMassErrorPowerLaw::nbodyHaloMassError_)
     select type (nbodyHaloMassError_)
     type is (nbodyHaloMassErrorPowerLaw)
        nbodyHaloMassError_=nbodyHaloMassErrorPowerLaw(                &
             &                                         stateVector(4), &
             &                                         stateVector(5)  &
             &                                        )
     end select
     haloMassFunctionConvolved_   =haloMassFunctionErrorConvolved(                              &
          &                                                       haloMassFunction_           , &
          &                                                       cosmologyParameters_        , &
          &                                                       nbodyHaloMassError_         , &
          &                                                       errorFractionalMaximum        &
          &                                                      )
     nullify(nbodyHaloMassError_)
     nullify(haloMassFunction_  )
     haloMassFunction_ => haloMassFunctionConvolved_
  case (haloMassFunctionErrorModelSphericalOverdensity)
     if (size(stateVector) /= 3 )                                                                  &
          & call Galacticus_Error_Report(                                                          & 
          &                              'likelihoodHaloMassFunctionEvaluate'                    , &
          &                              '3 parameters are required for this likelihood function'  &
          &                             )
     ! Use a mass function convolved with an error model for spherical overdensity algorithm errors.
     allocate(haloMassFunctionShethTormen::haloMassFunction_)
     select type (haloMassFunction_)
     type is (haloMassFunctionShethTormen)
        haloMassFunction_      =haloMassFunctionShethTormen(                              &
             &                                              cosmologyParameters_        , &
             &                                              cosmologicalMassVariance_   , &
             &                                              criticalOverdensity_        , &
             &                                              stateVector              (1), &
             &                                              stateVector              (2), &
             &                                              stateVector              (3)  &
             &                                             )
     end select
     darkMatterHaloScale_      => darkMatterHaloScale()
     darkMatterProfile_        => darkMatterProfile  ()
     allocate(nbodyHaloMassErrorSOHaloFinder::nbodyHaloMassError_)
     select type (nbodyHaloMassError_)
     type is (nbodyHaloMassErrorSOHaloFinder)
        nbodyHaloMassError_=nbodyHaloMassErrorSOHaloFinder(                      &
             &                                             darkMatterHaloScale_, &
             &                                             darkMatterProfile_  , &
             &                                             self%massParticle     &
             &                                            )
     end select
     haloMassFunctionConvolved_=haloMassFunctionErrorConvolved(                              &
          &                                                    haloMassFunction_           , &
          &                                                    cosmologyParameters_        , &
          &                                                    nbodyHaloMassError_         , &
          &                                                    errorFractionalMaximum        &
          &                                                   )
     nullify(darkMatterHaloScale_)
     nullify(darkMatterProfile_  )
     nullify(nbodyHaloMassError_ )
     nullify(haloMassFunction_   )
     haloMassFunction_ => haloMassFunctionConvolved_
  end select
  nullify(cosmologyParameters_     )
  nullify(cosmologicalMassVariance_)
  nullify(criticalOverdensity_     )
  ! Compute the mass function.
  allocate(massFunction(size(self%mass)))
  do i=1,size(self%mass)
     massFunction(i)=+haloMassFunction_%integrated(                     &
          &                                        self%time          , &
          &                                        self%massMinimum(i), &
          &                                        self%massMaximum(i)  &
          &                                       )                     &
          &          /log(                                              &
          &                                       +self%massMaximum(i)  &
          &                                       /self%massMinimum(i)  &
          &              )
  end do
  ! Evaluate the log-likelihood.
  difference                        =massFunction-self%massFunction
  likelihoodHaloMassFunctionEvaluate=-0.5d0*(difference*(self%inverseCovariance*difference))
  ! Clean up.
  deallocate(stateVector )
  deallocate(massFunction)
  return    
end function likelihoodHaloMassFunctionEvaluate

subroutine likelihoodHaloMassFunctionFunctionChanged(self)
  !% Respond to possible changes in the likelihood function.
  implicit none
  class(likelihoodHaloMassFunction), intent(inout) :: self
  !GCC$ attributes unused :: self

  return
end subroutine likelihoodHaloMassFunctionFunctionChanged

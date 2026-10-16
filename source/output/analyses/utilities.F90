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

!+    Contributions to this file made by: Claude.

!!{RST
Contains a module which provides a collection of utilities useful for on-the-fly analyses.
!!}

module Output_Analysis_Utilities
  !!{RST
  Provides a collection of utilities useful for on-the-fly analyses.
  !!}
  implicit none
  private
  public :: Output_Analysis_Output_Weight_Survey_Volume, Output_Analysis_Log_Likelihood_Normal, &
       &    Output_Analysis_Log_Likelihood_Relation    , Output_Analysis_Log_Likelihood_Write

contains

  function Output_Analysis_Output_Weight_Survey_Volume(surveyGeometry_,cosmologyFunctions_,outputTimes_,massLimit,magnitudeAbsoluteLimit,luminosity,allowSingleEpoch) result (outputWeight)
    !!{RST
    Compute output weights corresponding to the cosmological volumes associated with the given survey.
    !!}
    use            :: Cosmology_Functions, only : cosmologyFunctionsClass
    use            :: Display            , only : displayGreen           , displayReset
    use            :: Error              , only : Error_Report
    use            :: Geometry_Surveys   , only : surveyGeometryClass
    use, intrinsic :: ISO_C_Binding      , only : c_size_t
    use            :: ISO_Varying_String , only : operator(//)           , varying_string, assignment(=)
    use            :: Output_Times       , only : outputTimesClass
    implicit none
    double precision                         , dimension(:) , allocatable :: outputWeight
    class           (surveyGeometryClass    ), intent(inout)              :: surveyGeometry_
    class           (cosmologyFunctionsClass), intent(inout)              :: cosmologyFunctions_
    class           (outputTimesClass       ), intent(inout)              :: outputTimes_
    double precision                         , intent(in   ), optional    :: massLimit                   , magnitudeAbsoluteLimit, &
         &                                                                   luminosity
    logical                                  , intent(in   ), optional    :: allowSingleEpoch
    double precision                         , parameter                  :: timeTolerance        =1.0d-6
    integer         (c_size_t               )                             :: iOutput
    integer                                                               :: iField
    double precision                                                      :: timeMinimum                 , timeMaximum           , &
         &                                                                   distanceMinimum             , distanceMaximum       , &
         &                                                                   distanceSurveyMinimum       , distanceSurveyMaximum , &
         &                                                                   redshiftMinimum             , redshiftMaximum       , &
         &                                                                   time                        , timeMinimumFound      , &
         &                                                                   timeMaximumFound
    character       (len=12                 )                             :: redshiftLow                 , redshiftHigh
    type            (varying_string         )                             :: message
    !![
    <optionalArgument name="allowSingleEpoch" defaultsTo=".false." />
    !!]

    allocate(outputWeight(outputTimes_%count()))
    outputWeight=0.0d0
    if (outputTimes_%count() == 1_c_size_t) then
       ! Handle cases where we have just a single epoch output.
       if (allowSingleEpoch_) then
          ! Iterate over all fields.
          timeMinimumFound=+huge(0.0d0)
          timeMaximumFound=-huge(0.0d0)
          do iField=1,surveyGeometry_%fieldCount()
             ! Test whether the output epoch lies within the range of comoving distances for this field.
             time            =cosmologyFunctions_%cosmicTime            (cosmologyFunctions_%expansionFactorFromRedshift(outputTimes_%redshift(1_c_size_t)                                                         ))
             timeMinimum     =cosmologyFunctions_%timeAtDistanceComoving(surveyGeometry_    %distanceMaximum            (mass=massLimit,magnitudeAbsolute=magnitudeAbsoluteLimit,luminosity=luminosity,field=iField))
             timeMaximum     =cosmologyFunctions_%timeAtDistanceComoving(surveyGeometry_    %distanceMinimum            (mass=massLimit,magnitudeAbsolute=magnitudeAbsoluteLimit,luminosity=luminosity,field=iField))
             timeMinimumFound=min(timeMinimumFound,timeMinimum)
             timeMaximumFound=max(timeMaximumFound,timeMaximum)
             if     (                                           &
                  &   time >= timeMinimum*(1.0d0-timeTolerance) &
                  &  .and.                                      &
                  &   time <= timeMaximum*(1.0d0+timeTolerance) &
                  & ) outputWeight=1.0d0
          end do
          if (all(outputWeight == 0.0d0)) then
             message='zero weight for output times from survey geometry "'//surveyGeometry_%objectType()//'"'//char(10)
             write (redshiftLow ,'(f12.6)') cosmologyFunctions_%redshiftFromExpansionFactor(cosmologyFunctions_%expansionFactor(timeMaximumFound))
             write (redshiftHigh,'(f12.6)') cosmologyFunctions_%redshiftFromExpansionFactor(cosmologyFunctions_%expansionFactor(timeMinimumFound))
             message=message//'   required redshift: '//trim(redshiftLow)//' ≤ z ≤ '//trim(redshiftHigh)
             call Error_Report(message//{introspection:location})
          end if
       else
          timeMinimumFound=+huge(0.0d0)
          timeMaximumFound=-huge(0.0d0)
          do iField=1,surveyGeometry_%fieldCount()
             ! Test whether the output epoch lies within the range of comoving distances for this field.
             timeMinimum     =cosmologyFunctions_%timeAtDistanceComoving(surveyGeometry_%distanceMaximum(field=iField))
             timeMaximum     =cosmologyFunctions_%timeAtDistanceComoving(surveyGeometry_%distanceMinimum(field=iField))
             timeMinimumFound=min(timeMinimumFound,timeMinimum)
             timeMaximumFound=max(timeMaximumFound,timeMaximum)
          end do
          write (redshiftLow ,'(f12.6)') cosmologyFunctions_%redshiftFromExpansionFactor(cosmologyFunctions_%expansionFactor(timeMaximumFound))
          write (redshiftHigh,'(f12.6)') cosmologyFunctions_%redshiftFromExpansionFactor(cosmologyFunctions_%expansionFactor(timeMinimumFound))
          message='single epoch output not permitted'                                                                                     //char(10)// &
               &  displayGreen()//'HELP:'//displayReset()//' this output analysis requires at least two outputs in the redshift interval '          // &
               &  trim(adjustl(redshiftLow))                                                                                                        // &
               &  ' to '                                                                                                                            // &
               &  trim(adjustl(redshiftHigh))                                                                                                       // & 
               &  ' but fewer than two outputs exists in this range'                                                                      //char(10)// &
               &  '      add at least two outputs in this interval'                                                                       //char(10)// &
               &  '      this is required so that evolution over the redshift interval can be accounted for'
          call Error_Report(message//{introspection:location})
       end if
    else
       timeMinimumFound=+huge(0.0d0)
       timeMaximumFound=-huge(0.0d0)
       do iOutput=1,outputTimes_%count()
          do iField=1,surveyGeometry_%fieldCount()
             if (iOutput == outputTimes_%count()) then
                redshiftMinimum=     outputTimes_%redshift(iOutput)
             else
                redshiftMinimum=sqrt(outputTimes_%redshift(iOutput)*outputTimes_%redshift(iOutput+1))
             end if
             if (iOutput ==                              1) then
                redshiftMaximum=     outputTimes_%redshift(iOutput)
             else
                redshiftMaximum=sqrt(outputTimes_%redshift(iOutput)*outputTimes_%redshift(iOutput-1))
             end if
             timeMinimum          =    cosmologyFunctions_%cosmicTime      (cosmologyFunctions_%expansionFactorFromRedshift(redshiftMaximum))
             timeMaximum          =    cosmologyFunctions_%cosmicTime      (cosmologyFunctions_%expansionFactorFromRedshift(redshiftMinimum))
             distanceSurveyMinimum=    surveyGeometry_    %distanceMinimum (mass=massLimit  ,magnitudeAbsolute=magnitudeAbsoluteLimit,luminosity=luminosity,field=iField)
             distanceSurveyMaximum=    surveyGeometry_    %distanceMaximum (mass=massLimit  ,magnitudeAbsolute=magnitudeAbsoluteLimit,luminosity=luminosity,field=iField)
             distanceMinimum      =max(                                                                                                                                    &
                  &                    cosmologyFunctions_%distanceComoving(     timeMaximum                                                                            ), &
                  &                    distanceSurveyMinimum                                                                                                               &
                  &                   )
             distanceMaximum      =min(                                                                                                                                    &
                  &                    cosmologyFunctions_%distanceComoving(     timeMinimum                                                                            ), &
                  &                    distanceSurveyMaximum                                                                                                               &
                  &                   )
             timeMinimumFound     =min(timeMinimumFound,cosmologyFunctions_%timeAtDistanceComoving(distanceSurveyMaximum))
             timeMaximumFound     =max(timeMaximumFound,cosmologyFunctions_%timeAtDistanceComoving(distanceSurveyMinimum))
             outputWeight                      (iOutput)  &
                  & =outputWeight              (iOutput)  &
                  & +surveyGeometry_%solidAngle(iField )  &
                  & /3.0d0                                &
                  & *max(                                 &
                  &      +0.0d0                         , &
                  &      +distanceMaximum**3              &
                  &      -distanceMinimum**3              &
                  &    )
          end do
       end do
       where(outputWeight < 0.0d0)
          outputWeight=0.0d0
       end where
       if      (                                 &
            &    any(outputWeight >       0.0d0) &
            &  ) then
          outputWeight=+    outputWeight  &
               &       /sum(outputWeight)
       else if (                                 &
            &    timeMinimumFound < +huge(0.0d0) &
            &   .and.                            &
            &    timeMaximumFound > -huge(0.0d0) &
            &  ) then
          message='zero weight for all output times from survey geometry "'//surveyGeometry_%objectType()//'"'//char(10)
          write (redshiftLow ,'(f12.6)') cosmologyFunctions_%redshiftFromExpansionFactor(cosmologyFunctions_%expansionFactor(timeMaximumFound))
          write (redshiftHigh,'(f12.6)') cosmologyFunctions_%redshiftFromExpansionFactor(cosmologyFunctions_%expansionFactor(timeMinimumFound))
          message=message//'   required redshift(s): '//trim(redshiftLow)//' ≤ z ≤ '//trim(redshiftHigh)
          call Error_Report(message//{introspection:location})
       else
          message='zero weight for all output times from survey geometry "'//surveyGeometry_%objectType()//'"'
          call Error_Report(message//{introspection:location})
       end if
    end if
    return
  end function Output_Analysis_Output_Weight_Survey_Volume

  double precision function Output_Analysis_Log_Likelihood_Normal(difference,covarianceCombined,likelihoodNormalize) result(logLikelihood)
    !!{RST
    Compute the log-likelihood of a model given the difference, :math:`\Delta`, between model and target, and their combined
    covariance, :math:`C`, assuming a multivariate normal distribution:

    .. math::
     \log \mathcal{L} = -\frac{1}{2} \Delta^\mathrm{T} C^{-1} \Delta \left[ - \frac{1}{2} \log |C| - \frac{N}{2} \log(2\pi) \right],

    where the bracketed normalization terms are included only if ``likelihoodNormalize`` is true, and :math:`N` is the length
    of :math:`\Delta`. If the covariance can not be inverted the model is judged to be improbable.
    !!}
    use :: Interface_GSL               , only : GSL_Success
    use :: Linear_Algebra              , only : assignment(=), matrix, operator(*), vector
    use :: Models_Likelihoods_Constants, only : logImprobable
    use :: Numerical_Constants_Math    , only : Pi
    implicit none
    double precision        , dimension(:  ), intent(in   ) :: difference
    double precision        , dimension(:,:), intent(in   ) :: covarianceCombined
    logical                                 , intent(in   ) :: likelihoodNormalize
    type            (vector)                                :: residual
    type            (matrix)                                :: covariance
    integer                                                 :: status

    residual     =vector(difference        )
    covariance   =matrix(covarianceCombined)
    logLikelihood=-0.5d0*covariance%covarianceProduct(residual,status)
    if (status == GSL_Success) then
       if (likelihoodNormalize)                                        &
            & logLikelihood=+logLikelihood                             &
            &               -0.5d0*covariance%logarithmicDeterminant() &
            &               -0.5d0*dble(size(difference))*log(2.0d0*Pi)
    else
       logLikelihood       =+logImprobable
    end if
    return
  end function Output_Analysis_Log_Likelihood_Normal

  double precision function Output_Analysis_Log_Likelihood_Relation(value,covariance,valueTarget,covarianceTarget,likelihoodBins,likelihoodBinsAutomatic,valueMinimum,likelihoodNormalize) result(logLikelihood)
    !!{RST
    Compute the log-likelihood of a binned relation (e.g. the mean or scatter of some property as a function of another) given
    model and target values and covariances. The bins used are those listed in ``likelihoodBins`` (or all bins if that list is
    empty), unless ``likelihoodBinsAutomatic`` is true, in which case all bins with non-zero model value are used. If any
    active bin has a model value no greater than ``valueMinimum`` the model is judged to be improbable. Otherwise the
    log-likelihood is evaluated by ``Output_Analysis_Log_Likelihood_Normal``.
    !!}
    use, intrinsic :: ISO_C_Binding               , only : c_size_t
    use            :: Models_Likelihoods_Constants, only : logImprobable
    implicit none
    double precision          , dimension(:  ), intent(in   ) :: value                  , valueTarget
    double precision          , dimension(:,:), intent(in   ) :: covariance             , covarianceTarget
    integer         (c_size_t), dimension(:  ), intent(in   ) :: likelihoodBins
    logical                                   , intent(in   ) :: likelihoodBinsAutomatic, likelihoodNormalize
    double precision                          , intent(in   ) :: valueMinimum
    integer         (c_size_t), dimension(:  ), allocatable   :: bins
    double precision          , dimension(:,:), allocatable   :: covarianceCombined     , covarianceCombinedSelected
    double precision          , dimension(:  ), allocatable   :: difference             , differenceSelected
    integer                                                   :: i                      , j

    ! Determine which bins to use in the likelihood analysis.
    if (likelihoodBinsAutomatic) then
       j=0
       do i=1,size(value)
          if (value(i) /= 0.0d0) j=j+1
       end do
       allocate(bins(j))
       j=0
       do i=1,size(value)
          if (value(i) /= 0.0d0) then
             j=j+1
             bins(j)=i
          end if
       end do
    else
       allocate(bins,source=likelihoodBins)
    end if
    if     (                                                          &
         &   (size(bins) == 0 .and. any(value       <= valueMinimum)) &
         &  .or.                                                      &
         &                          any(value(bins) <= valueMinimum)  &
         & ) then
       ! If any active bins contain zero galaxies, judge this model to be improbable.
       logLikelihood=logImprobable
    else
       ! Compute difference with the target dataset.
       allocate(difference        ,mold=value     )
       allocate(covarianceCombined,mold=covariance)
       difference        =+value     -valueTarget
       covarianceCombined=+covariance+covarianceTarget
       ! Construct a reduced set of bins.
       if (size(bins) > 0) then
          allocate(differenceSelected        (size(bins)           ))
          allocate(covarianceCombinedSelected(size(bins),size(bins)))
          do i=1,size(bins)
             differenceSelected           (i  )=difference        (bins(i)        )
             do j=1,size(bins)
                covarianceCombinedSelected(i,j)=covarianceCombined(bins(i),bins(j))
             end do
          end do
       else
          allocate(differenceSelected        ,source=difference        )
          allocate(covarianceCombinedSelected,source=covarianceCombined)
       end if
       logLikelihood=Output_Analysis_Log_Likelihood_Normal(differenceSelected,covarianceCombinedSelected,likelihoodNormalize)
    end if
    return
  end function Output_Analysis_Log_Likelihood_Relation

  subroutine Output_Analysis_Log_Likelihood_Write(analysisLabel,logLikelihood,groupName,selection)
    !!{RST
    Write the log-likelihood of an analysis as an attribute of its group in the output file, overwriting any existing value. If
    present, ``selection`` is also written as an attribute of the group.
    !!}
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5Group
    use :: ISO_Varying_String, only : char      , varying_string
    use :: Output_HDF5       , only : outputFile
    implicit none
    type            (varying_string), intent(in   )           :: analysisLabel
    double precision                , intent(in   )           :: logLikelihood
    type            (varying_string), intent(in   ), optional :: groupName    , selection
    type            (hdf5Group     ), target                  :: analysesGroup, subGroup
    type            (hdf5Group     ), pointer                 :: inGroup
    type            (hdf5Group     )                          :: analysisGroup

    !$ call hdf5Access%set()
    analysesGroup =  outputFile   %openGroup('analyses'     )
    inGroup       => analysesGroup
    if (present(groupName)) then
       subGroup   =  analysesGroup%openGroup(char(groupName))
       inGroup    => subGroup
    end if
    analysisGroup=inGroup%openGroup(char(analysisLabel))
    call analysisGroup%writeAttribute(logLikelihood,'logLikelihood')
    if (present(selection)) call analysisGroup%writeAttribute(selection,'selection')
    !$ call hdf5Access%unset()
    return
  end subroutine Output_Analysis_Log_Likelihood_Write

end module Output_Analysis_Utilities
